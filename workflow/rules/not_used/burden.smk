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
    