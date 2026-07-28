rule collect_results:
    input:
      tree = outdir + '/phylogeny/snp_alignment.rooted.treefile',
      matrix5 = outdir + '/detettore/ALL_presence-absence.5prime.tsv',
      matrix3 = outdir + '/detettore/ALL_presence-absence.3prime.tsv',
      matrix_metadata  = outdir + '/detettore/ALL_presence-absence.metadata.tsv',
      refins = outdir + '/detettore/ALL_reference_insertions.tsv',
      anchormap = outdir + '/detettore/anchor_map.tsv',
      copy_numbers = outdir + '/detettore/ALL_copy_numbers.tsv',
    output: outdir + '/is6110.RData'
    params:
        metadata = config['metadata'],
        genemap = config['annotations']['mtbc0_h37rv_map'],
        dejesus = config['dejesus_orfs'],
        resistance = config['who_resistance']
    conda: '../envs/R.yml'
    script: 
        '../scripts/collect_results.R'


rule benchmarking:
    input: 
        benchmark_summary = outdir + '/benchmarking/detettore/benchmark_summary.tsv'
    output: outdir + '/1_benchmarking.html'
    params:
        assembly_metadata = config['benchmarking']['assembly_metadata']
    conda: '../envs/r_benchmarking.yml'
    script: '../notebooks/1_benchmarking.Rmd'


rule copy_numbers:
    input: 
        rdata = outdir + '/is6110.RData',
        treedata = outdir + '/phylogeny/treedata.CN_ASR.tsv'
    output:
        html = outdir + '/2_copy_numbers.html'
    conda: '../envs/R.yml'
    script: '../notebooks/2_copynumbers.Rmd'


rule birth_death:
    input:
        branchsummary_5 = outdir + '/detettore/AS_5prime.branch_summary.tsv',
        branchsummary_3 = outdir + '/detettore/AS_3prime.branch_summary.tsv',
        treedata = outdir + '/phylogeny/treedata.CN_ASR.tsv',
        rdata = outdir + '/is6110.RData'
    output: outdir + '/3_birth_death.html'
    conda: '../envs/R.yml'
    script: '../notebooks/3_birthdeath.Rmd'


rule niche_space:
    input: 
        rdata =outdir + '/is6110.RData',
        hotspot_summary = outdir + '/detettore/hotspot_summary.tsv'
    
    output: outdir + '/4_niche_space.html'
    params: 
        annot = config['annotations']['mtbc0'],
        dejesus = config['dejesus_orfs'],
        genemap = config['annotations']['mtbc0_h37rv_map']
    conda: '../envs/R.yml'
    script: '../notebooks/4_nichespace.Rmd'


rule niche_model:
    input: 
        subtree = outdir + '/niche_model/subsampled_tree.rooted.nex',
        abc_params = outdir + '/niche_model/abc_params.tsv',
        abc_summary = outdir + '/niche_model/abc_summaries.tsv',
        rdata = outdir + '/is6110.RData'
    output: outdir + '/5_niche_model.html'
    conda: '../envs/R.yml'
    script: '../notebooks/5_nichemodel.Rmd'