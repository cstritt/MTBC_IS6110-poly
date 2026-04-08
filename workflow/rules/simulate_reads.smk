rule simulate_reads_MSv3:
    input: 
        assembly = lambda wildcards: assemblies[wildcards.assembly],
    output: 
        r1 = config['outdir'] + '/art_reads/{assembly}.COV{coverage}.RL{readlength}.1.fq.gz',
        r2 = config['outdir'] + '/art_reads/{assembly}.COV{coverage}.RL{readlength}.2.fq.gz'
    params:
        assembly="{assembly}",
        coverage="{coverage}",
        readlength="{readlength}",
        outdir=config['outdir']
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
          -o {params.outdir}/art_reads/{params.assembly}.COV{params.coverage}.RL{params.readlength}.

        gzip {params.outdir}/art_reads/{params.assembly}.COV{params.coverage}.RL{params.readlength}.1.fq
        gzip {params.outdir}/art_reads/{params.assembly}.COV{params.coverage}.RL{params.readlength}.2.fq

        """


rule detettore_simreads:
    input: 
        r1 = config['outdir'] + '/art_reads/{assembly}.COV{coverage}.RL{readlength}.1.fq.gz',
        r2 = config['outdir'] + '/art_reads/{assembly}.COV{coverage}.RL{readlength}.2.fq.gz',
    output: 
        insertions = config['outdir'] + '/detettore/simulated/{assembly}.COV{coverage}.RL{readlength}.reference_insertions.tsv',
        anchors = config['outdir'] + '/detettore/simulated/{assembly}.COV{coverage}.RL{readlength}.anchors.tsv' 
    params:
        target = 'software/detettore6110/resources/is_targets/IS6110.fasta',
        pref = '{assembly}.COV{coverage}.RL{readlength}',
        outdir = directory(config['outdir'] + '/detettore/simulated')
    conda: '../envs/detettore6110.yml'
    shell:
        """
        python3 software/detettore6110/detettore6110.py {input.r1} {input.r2} \
         -t {params.target} \
         -o {params.outdir} \
         -p {params.pref}
        """
