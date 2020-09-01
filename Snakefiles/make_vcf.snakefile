
rule run_dnascope:
    input:
        bam = "Align/{id}.sort.rmdup.bam",
        ref = REF
    output:
        "Make_Vcf/step1_haplotyper/{id}_gatk.vcf"
    params:
        gatk_install = config['params']['gatk_install'],
        dbsnp = config['params']['dbsnp_path']
    benchmark:
        "Benchmarks/make_vcf.run_hpcaller.{id}.txt"
    run:
        command = ["{params.gatk_install}", "gatk", "HaplotypeCaller"
                   "-R", "{input.ref}", "-I", "{input.bam}"]
        if params.dbsnp:
            command.extend(["-D", "{params.dbsnp}"])
        shell(" ".join(command))


rule keep_pass_vars:
    input:
        "Make_Vcf/step1_haplotyper/{id}_gatk.vcf"
    output:
        "Make_Vcf/step1_haplotyper/{id}_gatk_pass_vars.vcf"
    benchmark:
        "Benchmarks/make_vcf.keep_pass_vars.{id}.txt"
    shell:
        """
        awk '($1~/^#/ || $7=="PASS" || $7=="."){{print}}' {input} > {output}
        """

rule select_snps:
    input:
        "Make_Vcf/step1_haplotyper/{id}_gatk.vcf"
    output:
        "Make_Vcf/step2_benchmarking/{id}.snp.vcf.gz"
    benchmark:
        "Benchmarks/make_vcf.select_snps.{id}.txt"
    shell:
        "/opt/cgi-tools/bin/bcftools view -O z --type snps {input} > {output}"


rule select_indels:
    input:
        "Make_Vcf/step1_haplotyper/{id}_gatk.vcf"
    output:
        "Make_Vcf/step2_benchmarking/{id}.indel.vcf.gz"
    benchmark:
        "Benchmarks/make_vcf.select_indels.{id}.txt"
    shell:
        "/opt/cgi-tools/bin/bcftools view -O z --type indels {input} > {output}"


rule index_snps:
    input:
        "Make_Vcf/step2_benchmarking/{id}.snp.vcf.gz"
    output:
        "Make_Vcf/step2_benchmarking/{id}.snp.vcf.gz.tbi"
    benchmark:
        "Benchmarks/make_vcf.index_snps.{id}.txt"
    shell:
        "/opt/cgi-tools/bin/tabix -p vcf -f {input}"


rule index_indels:
    input:
        "Make_Vcf/step2_benchmarking/{id}.indel.vcf.gz"
    output:
        "Make_Vcf/step2_benchmarking/{id}.indel.vcf.gz.tbi"
    benchmark:
        "Benchmarks/make_vcf.index_indels.{id}.txt"
    shell:
        "/opt/cgi-tools/bin/tabix -p vcf -f {input}"

rule eval_snps:
    input:
        snp = "Make_Vcf/step2_benchmarking/{}.snp.vcf.gz".format(config['samples']['id']),
        index = "Make_Vcf/step2_benchmarking/{}.snp.vcf.gz.tbi".format(config['samples']['id'])
    output:
        "Make_Vcf/step2_benchmarking/snp_compare/summary.txt"
    params:
        benchmark_snp = config['benchmark']['benchmark_snp'],
        ref_sdf = config['benchmark']['ref_sdf'],
        bedfile = config['benchmark']['bedfile'],
        direc = "Make_Vcf/step2_benchmarking/snp_compare"
    benchmark:
        "Benchmarks/make_vcf.eval_snps.txt"
    shell:
        "rm -r {params.direc};"
        "/research/rv-02/home/qmao/Scripts/rtg-tools-3.8.4/rtg vcfeval "
            "-b {params.benchmark_snp} "
            "-c {input.snp} "
            "-e {params.bedfile} "
            "-t {params.ref_sdf} "
            "-o {params.direc}"


rule eval_indels:
    input:
        indel = "Make_Vcf/step2_benchmarking/{}.indel.vcf.gz".format(config['samples']['id']),
        index = "Make_Vcf/step2_benchmarking/{}.indel.vcf.gz.tbi".format(config['samples']['id'])
    output:
        "Make_Vcf/step2_benchmarking/indel_compare/summary.txt"
    params:
        benchmark_indel = config['benchmark']['benchmark_indel'],
        ref_sdf = config['benchmark']['ref_sdf'],
        bedfile = config['benchmark']['bedfile'],
        direc = "Make_Vcf/step2_benchmarking/indel_compare"
    benchmark:
        "Benchmarks/make_vcf.eval_indels.txt"
    shell:
        "rm -r {params.direc};"
        "/research/rv-02/home/qmao/Scripts/rtg-tools-3.8.4/rtg vcfeval "
            "-b {params.benchmark_indel} "
            "-c {input.indel} "
            "-e {params.bedfile} "
            "-t {params.ref_sdf} "
            "-o {params.direc} "

