
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

