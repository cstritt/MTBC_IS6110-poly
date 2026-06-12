
rule create_gnr_file:
    input: config['gnumbers']
    output: 'config/datasets/gnrs.txt'
    shell:
        """ 
        cut -f1 {input} > {output}
        """

rule depth_files_tsv:
    input: config['gnumbers']
    output: outdir + '/alignment/depth_files.tsv'
    params:
        path_to_db = config['path_to_db']
    shell:
        """
        while read GNUM; do

            A="${{GNUM:0:3}}"
            B="${{GNUM:3:2}}"
            C="${{GNUM:5:2}}"
            DEPTH={params.path_to_db}/$A/$B/$C/$GNUM.depth.gz
            echo -e "${{GNUM}}\t${{DEPTH}}" >> {output}

        done < {input}
        """

rule vcf_txt:
    input: config['vcf']
    output: outdir + '/alignment/vcf_path.txt'
    shell:
        """
        echo {input} > {output}
        """

rule create_alignment:
    input:
        lva = 'software/large_variable_alignment/get_alignment.py',
        vcf = outdir + '/alignment/vcf_path.txt',
        depth = outdir + '/alignment/depth_files.tsv'
    output: 
        aln =  outdir + '/alignment/snp_alignment.fasta'
    params: 
        outdir = outdir + '/alignment',
        exclude = 'software/large_variable_alignment/resources/exclude/excluded_regions.bed'
    conda: '../envs/alignment.yml'
    threads: 20
    shell:
        """
        python3 {input.lva} \
          -i {input.vcf} \
          -d {input.depth} \
          -e {params.exclude} \
          -o {params.outdir} \
          -t {threads}
        """


rule rapidnj:
    input:  outdir + '/alignment/snp_alignment.fasta'
    output:  outdir + '/phylogeny/nj_tree.nex'
    conda: '../envs/rapidnj.yml'
    shell:
        """
        rapidnj {input} \
          -i fa \
          -a kim \
          -m 100000 \
          -c 20 \
          -b 100 \
          > {output}
        """

rule iqtree:
    input: 
        iqtree = 'software/iqtree-3.0.1-Linux/bin/iqtree3',
        aln = outdir + '/alignment/snp_alignment.fasta'
    output: outdir + '/phylogeny/snp_alignment.treefile'
    params:
        bootstraps = config['phylogeny']['bootstraps'],
        prefix = 'snp_alignment',
        outdir = outdir + '/phylogeny/'
    shell:
        """
        {input.iqtree} \
            -s  {input.aln} \
            -m GTR+ASC \
            --pathogen \
            --prefix {params.prefix} \
            -T AUTO \
            --threads-max 20 \
            --mem 400G

        mv {params.prefix}* {params.outdir}
        """

rule iqtree_maple:
    input: 
        iqtree = 'software/iqtree-3.0.1-Linux/bin/iqtree3',
        aln = outdir + '/alignment/snp_alignment.fasta'
    output: outdir + '/phylogeny/snp_alignment.maple.treefile'
    params:
        bootstraps = config['phylogeny']['bootstraps'],
        prefix = 'snp_alignment.maple',
        outdir = outdir + '/phylogeny/'
    shell:
        """
        {input.iqtree} \
            -s  {input.aln} \
            -m GTR+ASC \
            --pathogen \
            --alrt 1000 \
            --prefix {params.prefix} \
            -T AUTO \
            --threads-max 40 \
            --mem 200G \
            --safe

        mv {params.prefix}* {params.outdir}
        """

rule root_and_drop:
    input: 
        njtree = outdir + '/phylogeny/nj_tree.nex',
        iqtree = outdir + '/phylogeny/snp_alignment.treefile'
    output: 
        njtree = outdir + '/phylogeny/nj_tree.rooted.nex',
        iqtree = outdir + '/phylogeny/snp_alignment.rooted.treefile'

    params:
        outgroup = config['phylogeny']['outgroup']
    conda: '../envs/bio.yml'
    shell:
        """
        Rscript scripts/root_tree.r {input.njtree} {params.outgroup} {output.njtree}
        Rscript scripts/root_tree.r {input.iqtree} {params.outgroup} {output.iqtree}

        """

