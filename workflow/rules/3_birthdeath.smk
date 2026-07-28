rule ASR_birth_death:
    input:
        tree = outdir + '/phylogeny/snp_alignment.rooted.treefile',
        matrix5 = outdir + '/detettore/ALL_presence-absence.5prime.tsv',
        matrix3 = outdir + '/detettore/ALL_presence-absence.3prime.tsv'
    output: 
        as5prime = outdir + '/detettore/AS_5prime.branch_summary.tsv',
        as3prime = outdir + '/detettore/AS_3prime.branch_summary.tsv'
    params: 
        outpath = outdir + '/detettore'
    conda: '../envs/asr.yml'
    shell:
        """
        Rscript scripts/ancestral_state_reconstruction.births_and_deaths.R \
          {input.tree} \
          {input.matrix5} \
          {params.outpath} \
          AS_5prime

        Rscript scripts/ancestral_state_reconstruction.births_and_deaths.R \
          {input.tree} \
          {input.matrix3} \
          {params.outpath} \
          AS_3prime
        """



    