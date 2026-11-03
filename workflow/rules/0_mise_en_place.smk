configfile: 'config/config.yaml'

# Map assembly names to assembly paths
assemblies = {}
with open(config['benchmarking']['assemblies']) as f:
    next(f)
    for line in f:
        fields = line.strip().split("\t")
        assemblies[fields[0]] = fields[1]

# Create dictionary of G numbers with path to cram file as entries
gnrs = {}

# Parse RUN_ACCESSION from metadata for NCBI downloads
ncbi_accessions = {}
read_types = {}  # Store read type (PE vs SE) for each sample

with open(config["metadata"]) as f:
    header = next(f).strip().split("\t")
    gnumberi = header.index("GNUMBER")
    run_accessioni = header.index("RUN_ACCESSION")
    read_typei = header.index("READ TYPE") if "READ TYPE" in header else -1

    for line in f:
        fields = line.strip().split("\t")
        gnumber = fields[gnumberi]
        run_accession = fields[run_accessioni]
        if run_accession and run_accession != "":
            ncbi_accessions[gnumber] = run_accession
        if read_typei >= 0:
            read_types[gnumber] = fields[read_typei]

# Configure gnrs dictionary based on data_source
if config['data_source'] == "scicore":
    # Use existing CRAM paths from sciCORE
    with open(config['gnumbers']) as f:
        for line in f:
            gnr = line.strip().split('\t')[0]
            a = gnr[0:3]
            b = gnr[3:5]
            c = gnr[5:]
            gnrs[gnr] = f'{config['path_to_db']}/{a}/{b}/{c}/{gnr}.cram'
elif config['data_source'] == "ncbi":
    # Use downloaded FASTQ files from NCBI
    ncbi_path = config['ncbi_download_path']
    if not os.path.exists(ncbi_path):
        os.makedirs(ncbi_path)
    with open(config['gnumbers']) as f:
        for line in f:
            gnr = line.strip().split('\t')[0]
            a = gnr[0:3]
            b = gnr[3:5]
            c = gnr[5:]
            # Use first accession if multiple are separated by semicolon
            accession = ncbi_accessions.get(gnr, "").split(';')[0] if gnr in ncbi_accessions else ""
            if accession:
                # For FASTQ, we'll use the R1 file as the main input
                gnrs[gnr] = f'{ncbi_path}/{a}/{b}/{c}/{gnr}_1.fastq.gz'
            else:
                # Skip samples without NCBI accession
                continue
elif config['data_source'] == "db":
    # Extract FASTQ from IS6110 database
    import sqlite3
    db_path = config['db_path']
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT sample_name FROM samples ORDER BY sample_name;")
    for row in cursor.fetchall():
        gnr = row[0]
        a = gnr[0:3]
        b = gnr[3:5]
        c = gnr[5:]
        gnrs[gnr] = f'{outdir}/db_extracted/{a}/{b}/{c}/{gnr}.fastq.gz'
    conn.close()


# Output directory
import os
outdir = config['outdir']['main']
if not os.path.exists(outdir):
    os.makedirs(outdir)

outdir_bm = config['outdir']['benchmarking']
if not os.path.exists(outdir_bm):
    os.makedirs(outdir_bm)

# Create output folder structure for G numbers
detettore_gnr_subdirs= []
for gnr in gnrs.keys():
    detettore_gnr_subdirs.append(os.path.join(gnr[:3], gnr[3:5], gnr[5:]))

# Create lineage dictionary
lineages = {}

with open(config["metadata"]) as f:
    header = next(f).strip().split("\t")
    gnumberi = header.index("GNUMBER")
    lineagei = header.index("LINEAGE_x")

    for line in f:
        fields = line.strip().split("\t")
        gnumber = fields[gnumberi]
        lineage = fields[lineagei]
        if lineage not in lineages:
            lineages[lineage] = []
        lineages[lineage].append(gnumber)


rule clone_detettore:
    output: 'software/detettore6110/detettore6110.py'
    shell: 
        """
        git clone https://github.com/cstritt/detettore6110.git software/detettore6110
        touch {output}
        """

rule clone_lva:
    output: 
        'software/large_variable_alignment/get_alignment.py'
    shell: 
        """
        git clone https://github.com/cstritt/large_variable_alignment.git software/large_variable_alignment
        touch {output}
        
        """

rule clone_IS_simulation:
    output:
        'software/IS_simulation/abc_sampler.py'
    shell:
        """
        git clone https://github.com/cstritt/IS_simulation.git software/IS_simulation
        """

rule get_iqtree:
    output:
        'software/iqtree-3.0.1-Linux/bin/iqtree3'
    shell:
        """
        wget https://github.com/iqtree/iqtree3/releases/download/v3.0.1/iqtree-3.0.1-Linux.tar.gz software/

        tar -xzf software/iqtree-3.0.1-Linux.tar.gz
        rm software/iqtree-3.0.1-Linux.tar.gz

        """

# NCBI SRA download rules (only used when data_source == "ncbi")
if config['data_source'] == "ncbi":

    rule download_sra:
        input:
            lambda wildcards: ncbi_accessions[wildcards.gnr].split(';')[0] if wildcards.gnr in ncbi_accessions else ""
        output:
            sra = temp("{config[ncbi_download_path]}/{gsubdir}/{gnr}.sra")
        params:
            accession = lambda wildcards: ncbi_accessions[wildcards.gnr].split(';')[0] if wildcards.gnr in ncbi_accessions else "",
            outdir = lambda wildcards: f"{config['ncbi_download_path']}/{wildcards.gsubdir}"
        conda: '../envs/sratoolkit.yml'
        shell:
            """
            mkdir -p {params.outdir}
            prefetch {params.accession}
            mv $HOME/ncbi/public/sra/{params.accession}.sra {output.sra}
            """

    rule convert_to_fastq:
        input:
            sra = "{config[ncbi_download_path]}/{gsubdir}/{gnr}.sra"
        output:
            r1 = "{config[ncbi_download_path]}/{gsubdir}/{gnr}_1.fastq.gz",
            r2 = "{config[ncbi_download_path]}/{gsubdir}/{gnr}_2.fastq.gz"
        params:
            threads = config['sra_threads'],
            read_type = lambda wildcards: read_types.get(wildcards.gnr, "PE"),
            outdir = lambda wildcards: f"{config['ncbi_download_path']}/{wildcards.gsubdir}",
            gnr = "{gnr}"
        conda: '../envs/sratoolkit.yml'
        shell:
            """
            mkdir -p {params.outdir}
            if [ "{params.read_type}" == "PE" ]; then
                fasterq-dump --split-files --threads {params.threads} --gzip {input.sra} -O {params.outdir}/
                mv {params.outdir}/{params.gnr}_1.fastq.gz {output.r1}
                mv {params.outdir}/{params.gnr}_2.fastq.gz {output.r2}
            else
                fasterq-dump --threads {params.threads} --gzip {input.sra} -O {params.outdir}/
                mv {params.outdir}/{params.gnr}.fastq.gz {output.r1}
                touch {output.r2}
            fi
            """

    rule download_all_ncbi:
        input:
            expand("{config[ncbi_download_path]}/{gsubdir}/{gnr}_1.fastq.gz",
                   zip, gsubdir=detettore_gnr_subdirs, gnr=gnrs.keys())