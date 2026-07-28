rule subsample_tree:
    input: 
        tree = outdir + '/phylogeny/snp_alignment.rooted.treefile'
    output: outdir + '/niche_model/subsampled_tree.rooted.nex'
    params:
        metadata = config['metadata'],
        ntips = config['niche_model_abc']['ntips']
    conda: '../envs/R.yml'
    shell:
        """
        Rscript scripts/subsample_tree.R {input.tree} {params.metadata} {params.ntips} {output}
        """

rule run_abc:
    input: 
        tree = outdir + '/niche_model/subsampled_tree.rooted.nex',
    output: 
        parameters = outdir + '/niche_model/abc_params.tsv',
        summarystats = outdir + '/niche_model/abc_summaries.tsv'
    params: 
        metadata = config['metadata'],
        nsimulations = config['niche_model_abc']['nsimulations'],
        outpath = outdir + '/niche_model'
    conda: '../envs/niche_abc.yml'
    shell:
        """
        python scripts/niche_model_abc.py {input.tree} {params.nsimulations} {params.metadata} {params.outpath}
        """
