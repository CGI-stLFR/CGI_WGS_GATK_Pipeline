# Rule for mapping reads
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
        gatk_install = config['params']['gatk_install'],
        readgroup = r'@RG\tID:{0}\tSM:{0}\tPL:{1}'.format(config['samples']['id'],
                                                          config['params']['platform'])
    benchmark:
        "Benchmarks/main.map_reads.txt"
    shell:
        # map with readgroup and comments appended, this creates a BX tag with the barcode
        # Also includes sorting
        "bwa mem -M -R '{params.readgroup}' -C "
            "-t {threads} {input.ref} {input.fqs} 2>Align/aln.err | tee {output.sam} | "
        "samtools sort -o {output.bam} -@ {threads} -O bam -"


# Perform the deduplication step and generate metrics
rule mark_dups:
    input:
        bam = "Align/{id}.sort.bam",
    output:
        bam = "Align/{id}.sort.rmdup.bam",
        metrics = "Align/{id}_dedup_metrics.txt"
    params:
        gatk_install = config['params']['gatk_install']
    benchmark:
        "Benchmarks/main.mark_dups.{id}.txt"
    shell:
        "{params.gatk_install} MarkDuplicates -I {input.bam} "
            "-O {output.bam} "
            "-M {output.metrics}"


# Index markdups bam
rule index_mark_dups:
    input:
        bam = "Align/{id}.sort.rmdup.bam"
    output:
        bai = "Align/{id}.sort.rmdup.bam.bai"
    shell:
        "samtools index {input}"

# the mark duplicates rule just addds the flag, the below removes duplicates
# hapcut uses this sam file
# We also remove barcode uninformative reads (0_0_0)
rule remove_duplicates:
    input:
        "Align/{id}.sort.rmdup.bam"
    output:
        "Align/{id}.sort.removedup_rm000.sam"
    benchmark:
        "Benchmarks/calc_frag_len.remove_duplicates.{id}.txt"
    shell:
        "samtools view -h -F 0x400 {input} | "
        "awk -F $'\t' '($1!~/#0_0_0$/){{print}}' > {output}"

# This step parses the duplicate metrics and creates a more readable summary
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
