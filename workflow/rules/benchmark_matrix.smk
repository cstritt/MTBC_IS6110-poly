
rule minigraph_pangenome:
    input:
        assemblies = lambda wc: [assemblies[s] for s in sorted(assemblies)]
    output:
        gfa = outdir_bm + '/benchmarking/minigraph/pangenomegraph.gfa'
    threads: 8
    conda: '../envs/minigraph.yml'
    params:
        preset = "-cxggs"
    shell:
        """
        minigraph {params.preset} -t {threads} {input.assemblies} > {output.gfa}
        """

rule assembly_paths:
    input:
        assembly = lambda wc: assemblies[wc.assembly],
        gfa = outdir_bm + '/benchmarking/minigraph/pangenomegraph.gfa'
    output:
        bed = outdir_bm + '/benchmarking/minigraph/{assembly}.bed'
    threads: 8
    conda: '../envs/minigraph.yml'
    shell:
        """
        minigraph -cxasm --call -t {threads} {input.gfa} {input.assembly} > {output.bed}

        """

rule call_bubbles:
    input:
        gfa = outdir_bm + '/benchmarking/minigraph/pangenomegraph.gfa'
    output:
        bed = outdir_bm + '/benchmarking/minigraph/bubbles.bed'
    conda: '../envs/minigraph.yml'
    shell:
        """
        gfatools bubble {input.gfa} > {output.bed}
        """

rule find_IS_bubbles:
    input: 
        bubbles = outdir_bm + '/benchmarking/minigraph/bubbles.bed',
        is_seq = 'software/detettore6110/resources/is_targets/IS6110.fasta'
    output: outdir_bm + '/benchmarking/minigraph/IS6110_bubbles.tsv'
    conda: '../envs/minigraph.yml'
    shell:
        """
        scripts/find_IS_bubbles.py {input} {output}
        """





rule all:
    input:
        minigraph_pangenome = outdir_bm + '/benchmarking/minigraph/pangenomegraph.gfa',
        assembly_paths = outdir_bm + '/benchmarking/minigraph/{assembly}.bed',
        call_bubbles = outdir_bm + '/benchmarking/minigraph/bubbles.bed'
