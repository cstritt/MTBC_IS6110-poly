rule detettore:
    input: 
        reads = lambda wildcards: gnrs[wildcards.gnr],
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
        python3 software/detettore6110/detettore6110.py {input.reads} \
         -t {params.target} \
         -r {params.reference} \
         -a {params.annotation} \
         -o {params.outdir} \
         -p {params.gnr}
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
