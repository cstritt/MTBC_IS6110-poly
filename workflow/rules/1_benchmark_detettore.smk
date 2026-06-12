rule isescan:
    input: lambda wildcards: assemblies[wildcards.assembly]
    output: 
      d = directory(outdir + '/benchmarking/isescan/{assembly}'),
      gff = outdir + '/benchmarking/isescan/{assembly}/{assembly}.fna.gff',
      fna = outdir + '/benchmarking/isescan/{assembly}/{assembly}.fna.is.fna'
    params:
      assembly='{assembly}'
    conda: '../envs/isescan.yml'
    shell:
        """
        isescan.py \
          --seqfile {input} \
          --output {output.d} \
          --nthread 8

        mv {output.d}/{params.assembly}/* {output.d}/
        rm -rf {output.d}/{params.assembly} hmm proteome
        """

rule simulate_reads_MSv3:
    input:
        assembly = lambda wildcards: assemblies[wildcards.assembly]
    output:
        r1 = outdir_bm + '/benchmarking/art_reads/{assembly}/COV{coverage}.RL{readlength}.1.fq.gz',
        r2 = outdir_bm + '/benchmarking/art_reads/{assembly}/COV{coverage}.RL{readlength}.2.fq.gz'
    params:
        coverage='{coverage}',
        readlength='{readlength}',
        outdir=outdir_bm + '/benchmarking/art_reads/{assembly}'
    conda: '../envs/art.yml'
    shell:
        """
        art_illumina \
          -ss MSv3 \
          -i {input.assembly} \
          -p \
          -na \
          -l {params.readlength} \
          -f {params.coverage} \
          -m 400 \
          -s 20 \
          -o {params.outdir}/COV{params.coverage}.RL{params.readlength}.

        gzip {params.outdir}/COV{params.coverage}.RL{params.readlength}.1.fq
        gzip {params.outdir}/COV{params.coverage}.RL{params.readlength}.2.fq

        """

rule run_detettore:
    input:
        r1 = outdir_bm + '/benchmarking/art_reads/{assembly}/COV{coverage}.RL{readlength}.1.fq.gz',
        r2 = outdir_bm + '/benchmarking/art_reads/{assembly}/COV{coverage}.RL{readlength}.2.fq.gz'
    output:
        insertions = outdir_bm + '/benchmarking/detettore/{assembly}/{assembly}_COV{coverage}_RL{readlength}.reference_insertions.tsv',
        anchors = outdir_bm + '/benchmarking/detettore/{assembly}/{assembly}_COV{coverage}_RL{readlength}.anchors.tsv' 
    params:
        target = 'software/detettore6110/resources/is_targets/IS6110.fasta',
        reference = 'software/detettore6110/resources/reference/MTBC0_v1.1.fasta',
        annotation = 'software/detettore6110/resources/reference/MTBC0v1.1_PGAP_annot.gff',
        outdir = outdir_bm + '/benchmarking/detettore/{assembly}',
        pref = '{assembly}_COV{coverage}_RL{readlength}'

    conda: '../envs/detettore6110.yml'
    shell:
        """
        python3 software/detettore6110/detettore6110.py {input.r1} {input.r2} \
         -t {params.target} \
         -r {params.reference} \
         -a {params.annotation} \
         -o {params.outdir} \
         -p {params.pref}
        """


rule benchmark_detettore:
    input:
        assembly = lambda wildcards: assemblies[wildcards.assembly],
        isescan = outdir_bm + '/benchmarking/isescan/{assembly}/{assembly}.fna.is.fna',
        detettore = outdir_bm + '/benchmarking/detettore/{assembly}/{assembly}_COV{coverage}_RL{readlength}.anchors.tsv'
    output: outdir_bm + '/benchmarking/detettore/evaluation/{assembly}_COV{coverage}_RL{readlength}.eval.tsv'
    log: outdir_bm + '/benchmarking/detettore/evaluation/{assembly}_COV{coverage}_RL{readlength}.eval_detailed.tsv'
    params:
        target = 'software/detettore6110/resources/is_targets/IS6110.fasta',
        target_cluster = config['detettore']['is_cluster'],
        min_id = config['benchmarking']['min_id'],
        anchor_length = config['benchmarking']['anchor_length']
    conda: '../envs/detettore6110.yml'
    shell:
        """
        python3 scripts/benchmark.py \
          --detettore {input.detettore} \
          --isescan {input.isescan} \
          --assembly {input.assembly} \
          --min_id {params.min_id} \
          --anchor_length {params.anchor_length} \
          > {output}
        
        """

rule combine_results:
    input: expand(outdir_bm + '/benchmarking/detettore/evaluation/{assembly}_COV{coverage}_RL{readlength}.eval.tsv', assembly = assemblies.keys(), coverage = config['benchmarking']['coverages'], readlength = config['benchmarking']['readlengths'])
    output: outdir_bm + '/benchmarking/detettore/benchmark_summary.tsv'
    run:

        header = [
            'assembly', 'coverage', 'readlength', 
            'copy_number', 'paired_IRs', 'paired_IRs_wrong_seq', 'unpaired_IRL', 'unpaired_IRR',
            'side', 'true_positives', 'false_negatives', 'false_positives'
        ]
        
        with open(output[0], "w") as outhandle:

            outhandle.write('\t'.join(header)+'\n')

            for FILE in input:

                name_parsed = FILE.split('/')[-1].split('.')[0].split('_')
                assembly = name_parsed[0]
                coverage = name_parsed[1][3:]
                readlength = name_parsed[2][2:]

                with open(FILE) as f:
                    next(f)
                    for line in f:
                        row = line.strip().split('\t')
                        row = [assembly, coverage, readlength] + row
                        outhandle.write('\t'.join(row)+'\n')

        outhandle.close()


rule summarize:
    input: outdir_bm + '/benchmarking/detettore/benchmark_summary.tsv'
    output: outdir_bm + '/benchmarking/summary.html'
    conda: '../envs/R.yml'
    script: '../notebooks/1_benchmarking.Rmd'
