# at some point we were doing like three different fragment calculations and only using one
# I left them commented out just in case but we only ever looked at the one output by Sentieon's metrics

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


# rule for calculating metrics based on gatk's software
# it outputs gc bias, map quality distribution, quality distribution
# insert size and alignment metrics
rule calculate_metrics:
    input:
        bam = "Align/{id}.sort.bam",
        ref = config['params']['ref_fa']
    output:
        aln_sum = "Align/gatk_metrics_{id}.alignment_summary_metrics",
    threads:
        config['threads']['metrics']
    params:
        gatk_install = config['params']['gatk_install']
    benchmark:
        "Benchmarks/metrics.calculate_metrics.{id}.txt"
    shell:
        "{params.gatk_install} CollectMultipleMetrics "
            "-R {input.ref} "
            "-I {input.bam} "
            "-O Align/gatk_metrics_{wildcards.id}"


# This runs an analysis of duplicates
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


# This generates a plot of the duplicate analysis
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


# flagstat generates a summary based on sam flags
# This is used for the summary report
rule run_flagstat:
    input:
        "Align/{}.sort.rmdup.bam".format(config['samples']['id'])
    output:
        "Align/flagstat_metric.txt"
    benchmark:
        "Benchmarks/metrics.run_flagstat.txt"
    shell:
        "samtools flagstat {input} > {output}"


# this looks at things like mapping rate and duplicate rate, amongst others
# We also use it for the summary report
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


# This looks at the coverage across the genome, as well as percent coverage at particular depths (4X, 10X, 30X)
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


# This is one of the plots that isn't looked at frequently but can be turned on in the config file
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


# Wenlan had another way at looking at GC bias and this is it
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


# As we can have various outputs depending on the config settings
# this function alters the inputs to the summary report to make sure
# all the necessary files are generated
def summary_report_input(wildcards):
    summary_report_files = ["Align/coverage_depth.txt",
                            "Align/picard_align_metrics.txt",
                            "Align/gatk_metrics_{}.alignment_summary_metrics".format(config['samples']['id'])]

    # Add stLFR specific metrics that the summary report uses
    if config['modules']['stLFR']:
        calc_frag_file = ["Calc_Frag_Length/frag_length_distribution.pdf",
                          "Calc_Frag_Length/n_read_distribution.pdf",
                          "Calc_Frag_Length/frag_and_bc_summary.txt",
                          "Calc_Frag_Length/frags_per_bc.pdf"]


        for split_dist in config['calc_frag']['split_dist']:
            for outfile in calc_frag_file:
                parts = outfile.split("/")
                summary_report_files.append(parts[0] + "_" + str(split_dist) + "/" + parts[1])

    # Add phasing specific metrics that the summary report uses
    if config['modules']['phasing']:
        summary_report_files.append("Make_Vcf/step4_longhap/longhap_results.txt")
        summary_report_files.append("Make_Vcf/step3_hapcut/step4_compare_with_refphasing/hapcut_eval.txt")

    return summary_report_files


# rule to generate the summary report
# This report has pretty much everything that goes into the online spreadsheet
# It will search for all the potential files
# If it can't find a file it'll just emit a message to stderr and keep going
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
        "python3 {params.toolsdir}/tools/summary_report_v4.py | "
        "tee > {output}"
