# include the config file
configfile: "config.yaml"

# Define the ref based on the config file
# Sort of acts like a global variable so you don't need to always type the whole thing
REF = config['params']['ref_fa']

# This defines the logic for which chromosomes to use
# Check if the chroms are defined
if not config['samples']['chroms']:
    CHROMS = []
    # If they aren't defined try and scrape the fasta.fai index
    try:
        with open(REF+".fai", "r") as fasta_index:
            for line in fasta_index:
                CHROMS.append(line.strip().split()[0])

    except Exception as e:
        print(f"{REF}.fai cannot be opened for parsing.")
        print(f"Exception: {e}")
        sys.exit(1)
# if the chroms are defined just use those
else:
    CHROMS = config['samples']['chroms']


# define a function to return target files based on config settings
def run_all_input(wildcards):
    # These will always be returned and include the summary report, GC tables and graphs
    # Coverage, flagstat, and other metrics
    run_all_files = ["summary_report.txt",
                     "Align/GC_Table.csv",
                     "Align/GC_Distribution.png",
                     "Align/Normalized_GC_Coverage.png",
                     "Align/{}_dedup_metrics2.txt".format(config['samples']['id']),
                     "Align/flagstat_metric.txt",
                     "Align/gatk_metrics_{}.alignment_summary_metrics".format(config['samples']['id'])]

    # if stLFR is set to true add fragment calculation metrics
    if config['modules']['stLFR']:
        calc_frag_file = ["Calc_Frag_Length/frag_length_distribution.pdf",
                          "Calc_Frag_Length/n_read_distribution.pdf",
                          "Calc_Frag_Length/frag_and_bc_summary.txt",
                          "Calc_Frag_Length/frags_per_bc.pdf"]

        # add output directories for the various split distances
        for split_dist in config['calc_frag']['split_dist']:
            for outfile in calc_frag_file:
                parts = outfile.split("/")
                run_all_files.append(parts[0] + "_" + str(split_dist) + "/" + parts[1])

    # if variant_calling is true add VCF
    if config['modules']['variant_calling']:
        run_all_files.append("Make_Vcf/step1_haplotyper/{}_gatk.vcf".format(config['samples']['id']))

    # if benchmarking is set add benchmark summaries
    if config['modules']['benchmarking']:
        run_all_files.append("Make_Vcf/step2_benchmarking/snp_compare/summary.txt")
        run_all_files.append("Make_Vcf/step2_benchmarking/indel_compare/summary.txt")

    # if mate_pair_analysis is set add it's output
    if config['modules']['mate_pair_analysis']:
        run_all_files.append("Align/Distance_Between_Reads/mate_pairs_distance.pdf")

    # if read_overlap_analysis is set add it's output
    if config['modules']['read_overlap_analysis']:
        run_all_files.append("Align/Distance_Between_Reads/distance_between_reads.pdf")

    # if duplicate plot is set add it's output
    if config['modules']['duplicate_plot']:
        run_all_files.append("Benchmarks/metrics.duplicate_plot.txt")

    # if phasing is set add targets for hapcut and longhap analysis
    if config['modules']['phasing']:
        for CHR in CHROMS:
            run_all_files.append("Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{}_hapblock_{}".format(config['samples']['id'], CHR))
            run_all_files.append("Make_Vcf/step4_longhap/LibOutDir/{}.vcf".format(CHR))

        run_all_files.extend(["Make_Vcf/step4_longhap/longhap_results.txt",
                              "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/hapcut_eval.txt",
                              "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/{}_hapcut.phased.vcf".format(config['samples']['id']),
                              "Make_Vcf/step4_longhap/{}_gatk_pass_vars.vcf".format(config['samples']['id'])])

    return run_all_files


# rule run all, the files above are the targets for snakemake
rule run_all:
    input:
        run_all_input


# Include other modules and rules to run
# When possible the include statement points to the CGI_WGS_Pipeline
# This should limit redundancy when editing rules, hopefully
include: "/research/rv-02/home/eanderson/CGI_WGS_GATK_Pipeline/Snakefiles/calc_frag_len.smk"
include: "/research/rv-02/home/eanderson/CGI_WGS_Pipeline/Snakefiles/remove_sam_analysis.smk"
include: "/research/rv-02/home/eanderson/CGI_WGS_GATK_Pipeline/Snakefiles/make_vcf.smk"
include: "/research/rv-02/home/eanderson/CGI_WGS_GATK_Pipeline/Snakefiles/metrics.smk"
include: "/research/rv-02/home/eanderson/CGI_WGS_GATK_Pipeline/Snakefiles/stlfr.main.smk"
include: "/research/rv-02/home/eanderson/CGI_WGS_Pipeline/Snakefiles/splitreads.smk"
include: "/research/rv-02/home/eanderson/CGI_WGS_GATK_Pipeline/Snakefiles/phasing.smk"

# source the bash_profile setup for the pipeline
# shell.prefix('source /research/rv-02/home/eanderson/CGI_WGS_Pipeline/Data_and_Tools/bash_profile; ')

