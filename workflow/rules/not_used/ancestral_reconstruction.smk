
rule pastml_copy_number:
    input:
        tree = config['outdir'] + '/phylogeny/nj_tree.nex'
        table = config['outdir'] + '/detettore/copy_number.tsv'
    output:
        pass
    conda: '../envs/pastml.yml'
    shell:
        """
        pastml {input.tree} {input.table}

        """

rule pastml_insertion_sites:
    input:
        tree = config['outdir'] + '/phylogeny/nj_tree.nex'
        table = config['outdir'] + '/detettore/presence-absence.tsv'
    output:
        pass
    conda: '../envs/pastml.yml'
    shell:
        """
        pastml {input.tree} {input.table}

        """
