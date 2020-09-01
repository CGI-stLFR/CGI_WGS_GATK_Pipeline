rule map_reads:
    input:
        ref = REF,
        fqs = expand("{sample}", sample=config['samples']['fastq'])
    output:
        bam = "Align/{}.sort.bam".format(config['samples']['id']), # may want this to be a temp file
        sam = "Align/{}.aln_mem.sam".format(config['samples']['id'])
    threads:
        config['threads']['bwa']
    params:
        sen_install = config['params']['sentieon_install'],
        sen_license = config['params']['sentieon_license'],
        readgroup = r'@RG\tID:{0}\tSM:{0}\tPL:{1}'.format(config['samples']['id'],
                                                          config['params']['platform'])
    benchmark:
        "Benchmarks/main.map_reads.txt"
    shell:
        "{params.sen_install}/bin/bwa mem -M -R '{params.readgroup}' -C "
            "-t {threads} {input.ref} {input.fqs} 2>Align/aln.err | tee {output.sam} | "
        "{params.sen_install}/bin/sentieon util sort -o {output.bam} -t {threads} --sam2bam -i -"


rule locus_collector:
    input:
        "Align/{id}.sort.bam"
    output:
        "Align/{id}_score.txt"
    threads:
        config['threads']['bwa']
    params:
        sen_install = config['params']['sentieon_install']
    benchmark:
        "Benchmarks/main.locus_collector.{id}.txt"
    shell:
        "{params.sen_install}/bin/sentieon driver  -t {threads} "
            "-i {input} "
            "--algo LocusCollector "
            "--fun score_info {output}"


rule mark_dups:
    input:
        bam = "Align/{id}.sort.bam",
        score = "Align/{id}_score.txt"
    output:
        bam = "Align/{id}.sort.rmdup.bam",
        metrics = "Align/{id}_dedup_metrics.txt"
    threads:
        config['threads']['bwa']
    params:
        sen_install = config['params']['sentieon_install']
    benchmark:
        "Benchmarks/main.mark_dups.{id}.txt"
    shell:
        "{params.sen_install}/bin/sentieon driver  -t {threads} "
            "-i {input.bam} "
            "--algo Dedup "
            "--score_info {input.score} "
            "--metrics {output.metrics} {output.bam}"


rule mark_dups_txt:
    input:
        "Align/{id}_dedup_metrics.txt"
    output:
        "Align/{id}_dedup_metrics2.txt"
    params:
        toolsdir = config['params']['toolsdir']
    benchmark:
        "Benchmarks/main.mark_dups_txt.{id}.txt"
    shell:
        "perl {params.toolsdir}/tools/mark_dups_txt.pl {input} {output}"


