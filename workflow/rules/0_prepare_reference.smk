# Rules for preparing reference files
# These rules generate derived reference files from the raw MTBC0 and H37Rv inputs

rule mtbc0_to_h37rv_genemap:
    input:
        mtbc0_annot = config['annotations']['mtbc0'],
        mtbc0_genome = config['detettore']['reference'],
        h37rv_annot = config['annotations']['h37rv'],
        h37rv_genome = 'data/reference/H37Rv/GCF_000195955.2_ASM19595v2_genomic.fna'
    output:
        genemap = config['annotations']['mtbc0_h37rv_map']
    conda: '../envs/blast.yml'
    shell:
        """
        python3 workflow/scripts/compare_to_h37rv.py \
            --mtbc0_annot {input.mtbc0_annot} \
            --mtbc0_genome {input.mtbc0_genome} \
            --h37rv_annot {input.h37rv_annot} \
            --h37rv_genome {input.h37rv_genome} \
            --output {output.genemap}
        """

rule mtbc0_cds_intergenic_fasta:
    input:
        mtbc0_genome = config['detettore']['reference'],
        mtbc0_annot = config['annotations']['mtbc0']
    output:
        fasta = config['annotations']['mtbc0_region_seqs'],
        burden = 'data/reference/MTBC0/MTBC0.genic_n_intergenic_regions.burden.tsv'
    conda: '../envs/biopython.yml'
    shell:
        """
        python3 workflow/scripts/create_CDS_fasta.py \
            --input_fasta {input.mtbc0_genome} \
            --input_gff {input.mtbc0_annot} \
            --output_fasta {output.fasta} \
            --output_burden {output.burden}
        """

rule mtbc0_genomic_regions:
    input:
        mtbc0_annot = config['annotations']['mtbc0'],
        mtbc0_genome = config['detettore']['reference']
    output:
        regions = config['annotations']['mtbc0_regions']
    params:
        genome_length = 4411532
    conda: '../envs/biopython.yml'
    shell:
        """
        python3 workflow/scripts/parse_gff_regions.GC.py \
            --gff {input.mtbc0_annot} \
            --fasta {input.mtbc0_genome} \
            --out {output.regions} \
            --genome_length {params.genome_length}
        """
