rule subsample_tree:
    input: 
        tree = outdir + '/phylogeny/snp_alignment.rooted.treefile'
    output: outdir + '/niche_model/subsampled_tree.rooted.nex'
    params:
        metadata = config['metadata'],
        ntips = config['abc']['ntips']
    conda: '../envs/R.yml'
    shell:
        """
        Rscript scripts/subsample_tree.R {input.tree} {params.metadata} {params.ntips} {output}
        """

rule abc_basic:
    input: 
        tree = outdir + '/niche_model/subsampled_tree.rooted.nex',
        abc_sampler = 'software/IS_simulation/abc_sampler.py'
    output: 
        parameters = outdir + '/niche_model/basic/abc_params.tsv',
        summarystats = outdir + '/niche_model/basic/abc_summaries.tsv'
    params: 
        priors = config['abc']['priors'],
        metadata = config['metadata'],
        nsimulations = config['abc']['nsimulations'],
        outpath = outdir + '/niche_model/basic'
    conda: '../envs/niche_abc.yml'
    shell:
        """
        python {input.abc_sampler} \
          {input.tree} \
          {params.priors} \
          {params.metadata} \
          {params.outpath} \
          --nsim {params.nsimulations} \
          --model basic
        """

rule abc_target_pref:
    input: 
        tree = outdir + '/niche_model/subsampled_tree.rooted.nex',
        abc_sampler = 'software/IS_simulation/abc_sampler.py'
    output: 
        parameters = outdir + '/niche_model/tsp/abc_params.tsv',
        summarystats = outdir + '/niche_model/tsp/abc_summaries.tsv'
    params: 
        priors = config['abc']['priors'],
        metadata = config['metadata'],
        nsimulations = config['abc']['nsimulations'],
        outpath = outdir + '/niche_model/tsp'
    conda: '../envs/niche_abc.yml'
    shell:
        """
        python {input.abc_sampler} \
          {input.tree} \
          {params.priors} \
          {params.metadata} \
          {params.outpath} \
          --nsim {params.nsimulations} \
          --model target_site_prefs
          
        """

