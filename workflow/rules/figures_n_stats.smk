rule prologue:
    input:
      metadata = config['metadata']
      tree = outdir + '/phylogeny/snp_alignment.rooted.treefile',
      matrix5 = outdir + '/detettore/ALL_presence-absence.5prime.tsv',
      matrix3 = outdir + '/detettore/ALL_presence-absence.3prime.tsv',
      matrix_metadata  = outdir + '/detettore/ALL_presence-absence.metadata.tsv',
      refins = outdir + '/detettore/ALL_reference_insertions.tsv',
      anchormap = outdir + '/detettore/anchor_map.tsv',
      copy_numbers = outdir + '/detettore/ALL_copy_numbers.tsv'
    output: outdir + '/is6110.RData'
    conda: '../envs/R.yml'
    script: '../notebooks/0_prologue.Rmd'

rule copy_numbers:
    input: outdir + '/is6110.RData'
    output:
        html = outdir + '/2_copynumbers.html',
        treedata = outdir + '/phylogeny/treedata.CN_ASR.tsv'
    conda: '../envs/R.yml'
    script: '../notebooks/2_copynumbers.Rmd'

rule birth_death:
    input:
        rdata = outdir + '/is6110.RData'
        treedata = outdir + '/phylogeny/treedata.CN_ASR.tsv',
    output:
        html = outdir + '/3_birth_death.html'
    conda: '../envs/R.yml'
    script: '../notebooks/3_birthdeath.Rmd'

rule niche_space:
    input: outdir + '/is6110.RData'
    output: outdir + '/4_niche_space.html'
    conda: '../envs/R.yml'
    script: '../notebooks/4_nichespace.Rmd'

rule niche_model:
    input: outdir + '/is6110.RData'
    output: outdir + '/5_niche_model.html'
    conda: '../envs/R.yml'
    script: '../notebooks/5_nichemodel.Rmd'