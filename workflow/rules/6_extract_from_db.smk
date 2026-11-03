# Rules for extracting IS6110 partial hits from the database to FASTQ files
# These rules are only used when data_source == "db"

if config['data_source'] == "db":

    # Extract FASTQ for a single strain from the database
    rule extract_strain_fastq:
        input:
            db = config['db_path'],
            script = 'scripts/extract_db_to_fastq.py'
        output:
            fastq = outdir + '/db_extracted/{gsubdir}/{gnr}.fastq.gz'
        params:
            strain = '{gnr}'
        conda: '../envs/detettore6110.yml'
        shell:
            """
            python3 {input.script} {input.db} {params.strain} {output.fastq}
            """
