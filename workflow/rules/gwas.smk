
rule get_gene_coordinates:
    input: 
        mtbc0_annot = config['annotations']['mtbc0'],
        h37rv_annot = config['annotations']['h37rv'],
        annot_map = config['annotations']['mtbc0_h37rv_map']
    output: outdir + '/gwas/genes.burden_by_gene.tsv'
    run:
        import re
        go_terms = ['GO:0006310', 'GO:0006281', 'GO:0006260']

        # Grep GO terms in MTBC0 annotation
        r3_genes = []
        with open(input.mtbc0_annot) as f:
            for line in f:
                for term in go_terms:
                    if term in line and line not in r3_genes:
                        r3_genes.append(line)

        #%% Get corresponding orthologs in H37Rv
        genemap = {}
        with open(input.annot_map) as f:
            next(f)
            for line in f:
                fields = line.strip().split()
                genemap[fields[0]] = fields[1]

        pattern = r'ID=([^;]+)'

        rv_orthologs = []

        for line in r3_genes:
            match = re.search(pattern, line)
            gene_id = match.group(1)[4:]
            
            if gene_id in genemap:
                rv = genemap[gene_id]
                rv_orthologs.append(rv)
            else:  # Two integrases have no ortholog in H37Rv
                print(gene_id)

        #%% Write burden files for 1) R3 and 2) all genes
        pattern = r'locus_tag=([^;]+)'
        chromosome = 'MTB_anc'

        outhandle_all = open(output[0], 'w') 

        with open(input.h37rv_annot) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if fields[2] == 'gene':
                    match = re.search(pattern, fields[-1])
                    locus_tag = match.group(1)
                    start = fields[3]
                    end = fields[4]

                    outhandle_all.write(f'{locus_tag}\t{chromosome}:{start}-{end}\n')
                    
        outhandle_all.close()
                    

rule subset_vcf:
    input: config['vcf']
    output: outdir + '/gwas/non_synonymous.vcf.gz'
    params:
        tmp_vcf = outdir + '/gwas/non_synonymous.vcf'
    conda: '../envs/pyseer.yml'
    shell:
        """
        bcftools view {input} | grep -v "synonymous_variant" | grep -v "intergenic_region" > {params.tmp_vcf}
        bcftools view -O z {params.tmp_vcf} > {output}
        bcftools index {output}
        rm {params.tmp_vcf}
        """


rule get_phenotypes:
    input: outdir + '/detettore/ALL_copy_numbers.tsv'
    output: outdir + '/gwas/copy_number_phenotypes.tsv'
    shell:
        """
        cut -f1,3 {input} > {output}
        """

rule get_similarity_matrix:
    input: outdir + '/phylogeny/snp_alignment.rooted.treefile'
    output: outdir + '/gwas/iqtree.similarity.tsv'
    conda: '../envs/pyseer.yml'
    shell:
        """
        python scripts/phylogeny_distance.py {input} --lmm > {output}
        """  

rule prepare_IS_vcf:
    input: outdir + '/detettore/ALL_reference_insertions.vcf'
    output: outdir + '/gwas/ALL_reference_insertions.vcf.gz'
    conda: '../envs/pyseer.yml'
    shell:
        """
        bcftools view -O z {input} > {output}
        bcftools index {output}
        """  

rule burden_snpsindels:
    input:
        vcf = outdir + '/gwas/non_synonymous.vcf.gz',
        phenotypes = outdir + '/gwas/copy_number_phenotypes.tsv',
        similarity = outdir + '/gwas/iqtree.similarity.tsv',
        burden_regions = outdir + '/gwas/genes.burden_by_gene.tsv'
    output: outdir + '/gwas/burden_snps_indels.tsv'
    conda: '../envs/pyseer.yml'

    shell:
        """
        pyseer \
            --lmm --print-samples \
            --cpu 20 \
            --vcf {input.vcf} \
            --phenotypes {input.phenotypes} \
            --similarity {input.similarity} \
            --burden {input.burden_regions} > {output}
        
        """


rule burden_is6110:
    input: 
        vcf = outdir + '/gwas/ALL_reference_insertions.vcf.gz',
        phenotypes = outdir + '/gwas/copy_number_phenotypes.tsv',
        similarity = outdir + '/gwas/iqtree.similarity.tsv'
    output: outdir + '/gwas/burden_is6110.tsv'
    params: 
        burden_regions = config['annotations']['mtbc0_regions']
    conda: '../envs/pyseer.yml'
    shell:
        """
        # Burden IS6110
        pyseer \
            --lmm --print-samples \
            --cpu 20 \
            --vcf {input.vcf} \
            --phenotypes {input.phenotypes} \
            --similarity {input.similarity} \
            --burden {params.burden_regions} > {output}
        """

rule permute_phenotypes:
    input: outdir + '/gwas/copy_number_phenotypes.tsv'
    output: outdir + '/gwas/permutations/permuted_data{i}.tsv'
    params:
        i = '{i}',
        outpath = outdir + '/gwas/permutations'
    run:
        import random

        strains = []
        copy_numbers = []

        i = params.i

        with open(input[0]) as f:
            header = next(f)
            for line in f:
                strain, cn = line.strip().split()
                strains.append(strain)
                copy_numbers.append(cn)

        # Permute copy numbers and write to file
        outhandle = open(f'{params.outpath}/permuted_data{i}.tsv', 'w')
        random.shuffle(copy_numbers)
        outhandle.write(header)
        for j in range(len(copy_numbers)):
            outhandle.write(f'{strains[j]}\t{copy_numbers[j]}\n')
        outhandle.close()
        

rule permutations_snpsindels:
    input:
        vcf = outdir + '/gwas/non_synonymous.vcf.gz',
        phenotypes = outdir + '/gwas/permutations/permuted_data{i}.tsv',
        similarity = outdir + '/gwas/iqtree.similarity.tsv',
        burden_all = outdir + '/gwas/genes.burden_by_gene.tsv'
    output: 
        burden_all = outdir + '/gwas/permutations/burden_snps_indels.{i}.tsv'
    conda: '../envs/pyseer.yml'
    shell:
        """
        pyseer \
            --lmm \
            --cpu 40 \
            --vcf {input.vcf} \
            --phenotypes {input.phenotypes} \
            --similarity {input.similarity} \
            --burden {input.burden_all} > {output.burden_all}
        """

rule permutations_is6110:
    input:
        vcf = outdir + '/gwas/ALL_reference_insertions.vcf.gz',
        phenotypes = outdir + '/gwas/permutations/permuted_data{i}.tsv',
        similarity = outdir + '/gwas/iqtree.similarity.tsv',
    output: 
        burden_all = outdir + '/gwas/permutations/burden_is6110.{i}.tsv'
    params: 
        burden_regions = config['annotations']['mtbc0_regions']
    conda: '../envs/pyseer.yml'
    shell:
        """
        pyseer \
            --lmm \
            --cpu 40 \
            --vcf {input.vcf} \
            --phenotypes {input.phenotypes} \
            --similarity {input.similarity} \
            --burden {params.burden_regions} > {output.burden_all}
        """