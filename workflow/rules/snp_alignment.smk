rule create_gnr_file:
    input: config['gnumbers']
    output: 'config/datasets/gnrs.txt'
    shell:
        """ 
        cut -f1 {input} > {output}
        """

rule depth_files_tsv:
    input: config['gnumbers']
    output: outdir + '/alignment/depth_files.tsv'
    params:
        path_to_db = config['path_to_db']
    shell:
        """
        while read GNUM; do

            A="${{GNUM:0:3}}"
            B="${{GNUM:3:2}}"
            C="${{GNUM:5:2}}"
            DEPTH={params.path_to_db}/$A/$B/$C/$GNUM.depth.gz
            echo -e "${{GNUM}}\t${{DEPTH}}" >> {output}

        done < {input}
        """

rule vcf_txt:
    input: config['vcf']
    output: outdir + '/alignment/vcf_path.txt'
    shell:
        """
        echo {input} > {output}
        """

rule create_alignment:
    input:
        lva = 'software/large_variable_alignment/get_alignment.py',
        vcf = outdir + '/alignment/vcf_path.txt',
        depth = outdir + '/alignment/depth_files.tsv'
    output: 
        aln =  outdir + '/alignment/snp_alignment.fasta'
    params: 
        outdir = outdir + '/alignment',
        exclude = 'software/large_variable_alignment/resources/exclude/excluded_regions.bed'
    conda: '../envs/alignment.yml'
    threads: 20
    shell:
        """
        python3 {input.lva} \
          -i {input.vcf} \
          -d {input.depth} \
          -e {params.exclude} \
          -o {params.outdir} \
          -t {threads}
        """
