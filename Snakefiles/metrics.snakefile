
#rule gc_depth:
#    input:
#        bam = "Align/{}.sort.rmdup.bam".format(config['samples']['id']),
#        fasta_gc = config['params']['gc_bias_index']
#    output:
#        exome = "Align/exome_gc_depth.txt",
#        genome = "Align/genome_gc_depth.txt"
#    params:
#        toolsdir = config['params']['toolsdir']
#    shell:
#        "samtools depth {input.bam} | "
#        "{params.toolsdir}/tools/gcBiasBinCaler cal {input.fasta_gc} - 500 "
#            "{output.exome} > {output.genome} 2> /dev/null"
#
#
#rule genome_gc_bias:
#    input:
#        "Align/genome_gc_depth.txt"
#    output:
#        "Align/genome_gc_bias.txt"
#    params:
#        toolsdir = config['params']['toolsdir'],
#        python = config['params']['wenlan_python']
#    shell:
#        "{params.python} {params.toolsdir}/tools/calGcBiasCurve.py {input} > {output}"
#
#
#rule exome_gc_bias:
#    input:
#        "Align/exome_gc_depth.txt"
#    output:
#        "Align/exome_gc_bias.txt"
#    params:
#        toolsdir = config['params']['toolsdir'],
#        python = config['params']['wenlan_python']
#    shell:
#        "{params.python} {params.toolsdir}/tools/calGcBiasCurve.py {input} > {output}"
#
#
#rule plot_genome_gc_bias:
#    input:
#        "Align/genome_gc_bias.txt"
#    output:
#        "Align/genome_gc_bias.pdf"
#    params:
#        toolsdir = config['params']['toolsdir']
#    shell:
#        "{params.toolsdir}/tools/R-3.2.1/bin/Rscript "
#        "{params.toolsdir}/tools/AlignmentGCbias.R "
#        "{input} {output} 2"
#
#
#rule plot_exome_gc_bias:
#    input:
#        "Align/exome_gc_bias.txt"
#    output:
#        "Align/exome_gc_bias.pdf"
#    params:
#        toolsdir = config['params']['toolsdir']
#    shell:
#        "{params.toolsdir}/tools/R-3.2.1/bin/Rscript "
#        "{params.toolsdir}/tools/AlignmentGCbias.R "
#        "{input} {output} 2"


rule calculate_metrics:
    input:
        bam = "Align/{id}.sort.bam",
        ref = config['params']['ref_fa']
    output:
        gc_sum = "Align/sentieon_gc_{id}_sum.txt",
        gc_met = "Align/sentieon_gc_{id}_metric.txt",
        mq = "Align/sentieon_mq_{id}_metric.txt",
        qd = "Align/sentieon_qd_{id}_metric.txt",
        insert = "Align/sentieon_is_{id}_metric.txt",
        aln = "Align/sentieon_aln_{id}_metric.txt"
    threads:
        config['threads']['metrics']
    params:
        sen_install = config['params']['sentieon_install']
    benchmark:
        "Benchmarks/metrics.calculate_metrics.{id}.txt"
    shell:
        "{params.sen_install}/bin/sentieon driver "
            "-t {threads} "
            "-r {input.ref} "
            "-i {input.bam} "
            "--algo GCBias --summary {output.gc_sum} {output.gc_met} "
            "--algo MeanQualityByCycle {output.mq} "
            "--algo QualDistribution {output.qd} "
            "--algo InsertSizeMetricAlgo {output.insert} "
            "--algo AlignmentStat {output.aln}"


rule plot_metrics:
    input:
        gc_met = "Align/sentieon_gc_{id}_metric.txt",
        mq = "Align/sentieon_mq_{id}_metric.txt",
        qd = "Align/sentieon_qd_{id}_metric.txt",
        insert = "Align/sentieon_is_{id}_metric.txt"
    output:
        "Align/sentieon_metrics_{id}.pdf"
    params:
        sen_install = config['params']['sentieon_install']
    benchmark:
        "Benchmarks/metrics.plot_metrics.{id}.txt"
    shell:
        "{params.sen_install}/bin/sentieon plot metrics "
            "-o {output} "
            "gc={input.gc_met} "
            "mq={input.mq} "
            "qd={input.qd} "
            "isize={input.insert}"


rule duplicate_analysis:
    input:
        "Align/{}.sort.rmdup.bam".format(config['samples']['id'])
    output:
        "Align/Duplicate_Analysis/dup_info",
        "Align/Duplicate_Analysis/duplicate_rate",
        "Align/Duplicate_Analysis/PE_dup_reads"
    params:
        toolsdir = config['params']['toolsdir']
    benchmark:
        "Benchmarks/metrics.duplicate_analysis.txt"
    shell:
        "samtools view {input} | "
            "perl {params.toolsdir}/tools/Duplicate_analysis.pl - Align/Duplicate_Analysis"


rule duplicate_plot:
    input:
        "Align/Duplicate_Analysis/dup_info"
    output:
        "Align/Duplicate_Analysis/duplicate.pdf"
    params:
        toolsdir = config['params']['toolsdir']
    benchmark:
        "Benchmarks/metrics.duplicate_plot.txt"
    shell:
        "{params.toolsdir}/tools/duplicate_statistics_new.R {input} {output}"


rule run_flagstat:
    input:
        "Align/{}.sort.rmdup.bam".format(config['samples']['id'])
    output:
        "Align/flagstat_metric.txt"
    benchmark:
        "Benchmarks/metrics.run_flagstat.txt"
    shell:
        "samtools flagstat {input} > {output}"


rule picard_align_metrics:
    input:
        "Align/{}.sort.rmdup.bam".format(config['samples']['id'])
    output:
        "Align/picard_align_metrics.txt"
    params:
        toolsdir = config['params']['toolsdir']
    benchmark:
        "Benchmarks/metrics.picard_align_metrics.txt"
    shell:
        "perl {params.toolsdir}/tools/picard.pl {input} "
        "{params.toolsdir}/tools/samtools-0.1.18/samtools > {output}"


rule coverage_depth:
    input:
        "Align/{}.sort.rmdup.bam".format(config['samples']['id'])
    output:
        "Align/coverage_depth.txt"
    params:
        toolsdir = config['params']['toolsdir'],
        ref = config['params']['ref_fa']
    benchmark:
        "Benchmarks/metrics.coverage_depth.txt"
    shell:
        "perl {params.toolsdir}/tools/depthV2.0.pl -l $({params.toolsdir}/tools/fasta_non_gapped_bases.py {params.ref}) {input} Align > {output}"


rule coverage_plot_sam:
    input:
        "Calc_Frag_Length/step1_removedup_rm000/{id}.sort.removedup_rm000.sam"
    output:
        "Align/Coverage_Plot/Sam_File/{id}_aa.sam"
    params:
        id = config['samples']['id']
    benchmark:
        "Benchmarks/metrics.coverage_plot_sam.{id}.txt"
    shell:
        "mkdir -p Align/Coverage_Plot/Sam_File; cd Align/Coverage_Plot/Sam_File; "
        "split -l 1000000 --additional-suffix=.sam ../../../{input} {params.id}_"


rule moar_gc_plots:
    input:
        sam = "Align/{}.sort.removedup_rm000.sam".format(config['samples']['id']),
        ref = "{}/data/Human_GC_Bins.csv".format(config['params']['toolsdir'])
    output:
        "Align/GC_Table.csv",
        "Align/GC_Distribution.png",
        "Align/Normalized_GC_Coverage.png"
    params:
        python = config['params']['gcbias_python'],
        toolsdir = config['params']['toolsdir']
    benchmark:
        "Benchmarks/metrics.moar_gc_plots.txt"
    shell:
        "{params.python} {params.toolsdir}/tools/GC_bias_20181217.py "
            "{input.sam} "
            "-r {input.ref} "
            "-o Align/"


def summary_report_input(wildcards):
    summary_report_files = ["Align/coverage_depth.txt",
                            "Align/picard_align_metrics.txt",
                            "Align/sentieon_is_{}_metric.txt".format(config['samples']['id'])]

    if config['modules']['stLFR']:
        calc_frag_file = ["Calc_Frag_Length/frag_length_distribution.pdf",
                          "Calc_Frag_Length/n_read_distribution.pdf",
                          "Calc_Frag_Length/frag_and_bc_summary.txt",
                          "Calc_Frag_Length/frags_per_bc.pdf"]


        for split_dist in config['calc_frag']['split_dist']:
            for outfile in calc_frag_file:
                parts = outfile.split("/")
                summary_report_files.append(parts[0] + "_" + str(split_dist) + "/" + parts[1])

    if config['modules']['phasing']:
        summary_report_files.append("Make_Vcf/step4_longhap/longhap_results.txt")
        summary_report_files.append("Make_Vcf/step3_hapcut/step4_compare_with_refphasing/hapcut_eval.txt") 

    return summary_report_files


rule generate_summary_report:
    input:
        summary_report_input
    output:
        "summary_report.txt"
    params:
        toolsdir = config['params']['toolsdir'],
        samp = config['samples']['id'],
        read_length = config['params']['read_len'],
        min_frag = config['calc_frag']['min_frag']
    benchmark:
        "Benchmarks/metrics.generate_summary_report.txt"
    shell:
        "python {params.toolsdir}/tools/summary_report_v3.py {params.samp} {params.read_length} {params.min_frag}| "
        "tee > {output}"

