import os

# Parse read types from metadata for NCBI mode
read_types_detettore = {}
if config['data_source'] == "ncbi":
    with open(config["metadata"]) as f:
        header = next(f).strip().split("\t")
        gnumberi = header.index("GNUMBER")
        read_typei = header.index("READ TYPE") if "READ TYPE" in header else -1

        for line in f:
            fields = line.strip().split("\t")
            gnumber = fields[gnumberi]
            if read_typei >= 0:
                read_types_detettore[gnumber] = fields[read_typei]

def get_reads_input(wildcards):
    """Return appropriate read input files based on data source and read type"""
    reads = gnrs[wildcards.gnr]
    if config['data_source'] == "ncbi":
        # For NCBI, check if paired-end
        is_paired = read_types_detettore.get(wildcards.gnr, "PE") == "PE"
        if is_paired:
            reads2 = reads.replace('_1.fastq.gz', '_2.fastq.gz')
            return reads, reads2
    # For scicore or single-end, return empty string for reads2
    return reads, []

# Define detettore rule
rule detettore:
    input:
        reads = lambda wildcards: get_reads_input(wildcards)[0],
        reads2 = lambda wildcards: get_reads_input(wildcards)[1],
        detettore = 'software/detettore6110/detettore6110.py'
    output:
        insertions = outdir + '/detettore/{gsubdir}/{gnr}.reference_insertions.tsv',
        anchors = outdir + '/detettore/{gsubdir}/{gnr}.anchors.tsv'
    params:
        target = 'software/detettore6110/resources/is_targets/IS6110.fasta',
        reference = 'software/detettore6110/resources/reference/MTBC0_v1.1.fasta',
        annotation = 'software/detettore6110/resources/reference/MTBC0v1.1_PGAP_annot.gff',
        gnr = '{gnr}',
        outdir = directory(outdir + '/detettore/{gsubdir}')
    conda: '../envs/detettore6110.yml'
    shell:
        """
        if [ -n "{input.reads2}" ] && [ -f "{input.reads2}" ]; then
            python3 software/detettore6110/detettore6110.py {input.reads} {input.reads2} \
             -t {params.target} \
             -r {params.reference} \
             -a {params.annotation} \
             -o {params.outdir} \
             -p {params.gnr}
        else
            python3 software/detettore6110/detettore6110.py {input.reads} \
             -t {params.target} \
             -r {params.reference} \
             -a {params.annotation} \
             -o {params.outdir} \
             -p {params.gnr}
        fi
        """

rule detettore_summarize_results:
    input: 
        anchors = expand(outdir + '/detettore/{gsubdir}/{gnr}.anchors.tsv', zip, gsubdir=detettore_gnr_subdirs, gnr=gnrs.keys()),
        insertions = expand(outdir + '/detettore/{gsubdir}/{gnr}.reference_insertions.tsv', zip, gsubdir=detettore_gnr_subdirs, gnr=gnrs.keys()),
        samples = 'config/datasets/gnrs.txt',
        detettore = 'software/detettore6110/detettore6110.py'
    output: 
        anchors = outdir + '/detettore/ALL_anchors.tsv',
        refins = outdir + '/detettore/ALL_reference_insertions.tsv',
        refins_vcf = outdir + '/detettore/ALL_reference_insertions.vcf',
        matrix5 = outdir + '/detettore/ALL_presence-absence.5prime.tsv',
        matrix3 = outdir + '/detettore/ALL_presence-absence.3prime.tsv',
        copynumbers = outdir + '/detettore/ALL_copy_numbers.tsv',
        metadata = outdir + '/detettore/ALL_presence-absence.metadata.tsv',
        ua_fasta = outdir + '/detettore/ALL_unique_anchors.fasta',
        anchor_map = outdir + '/detettore/anchor_map.tsv'
    params:
        outdir = outdir + '/detettore',
        seed_len = config['detettore']['seed_len'],
        maxmissing = config['detettore']['maxmissing'],
        reference = 'software/detettore6110/resources/reference/MTBC0_v1.1.fasta'


    conda: '../envs/detettore6110.yml'
    shell:
        """
        python3 software/detettore6110/summarize.py \
          -i {params.outdir} \
          -x {input.samples} \
          -r {params.reference} \
          -o {params.outdir} \
          -s {params.seed_len} \
          -mm {params.maxmissing}

        """


rule summarize_hotspots:
    input: 
        tree = outdir + '/phylogeny/snp_alignment.rooted.treefile',
        matrix5 = outdir + '/detettore/ALL_presence-absence.5prime.tsv',
        matrix_meta = outdir + '/detettore/ALL_presence-absence.metadata.tsv'
    output:
        hotspot_summary = outdir + '/detettore/hotspot_summary.tsv'
    conda: '../envs/asr.yml'
    shell:
        """
        Rscript scripts/infer_hotspots.R \
          {input.tree} \
          {input.matrix5} \
          {input.matrix_meta} \
          {output.hotspot_summary}
        """
