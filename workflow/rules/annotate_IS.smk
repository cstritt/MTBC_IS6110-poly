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
      