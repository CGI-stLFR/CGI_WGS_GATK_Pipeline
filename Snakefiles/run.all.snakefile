configfile: "config.yaml"


REF = config['params']['ref_fa']
if not config['samples']['chroms']:
    CHROMS = []
    try:
        with open(REF+".fai", "r") as fasta_index:
            for line in fasta_index:
                CHROMS.append(line.strip().split()[0])

    except Exception as e:
        print(f"{REF}.fai cannot be opened for parsing.")
        print(f"Exception: {e}")
        sys.exit(1)
else:
    CHROMS = config['samples']['chroms']

def run_all_input(wildcards):
    run_all_files = ["summary_report.txt",
                     "Align/GC_Table.csv",
                     "Align/GC_Distribution.png",
                     "Align/Normalized_GC_Coverage.png",
                     "Align/{}_dedup_metrics2.txt".format(config['samples']['id']),
                     "Align/flagstat_metric.txt",
                     "Align/sentieon_metrics_{}.pdf".format(config['samples']['id'])]

    if config['modules']['stLFR']:
        calc_frag_file = ["Calc_Frag_Length/frag_length_distribution.pdf",
                          "Calc_Frag_Length/n_read_distribution.pdf",
                          "Calc_Frag_Length/frag_and_bc_summary.txt",
                          "Calc_Frag_Length/frags_per_bc.pdf"]

        for split_dist in config['calc_frag']['split_dist']:
            for outfile in calc_frag_file:
                parts = outfile.split("/")
                run_all_files.append(parts[0] + "_" + str(split_dist) + "/" + parts[1])

    if config['modules']['variant_calling']: 
        run_all_files.append("Make_Vcf/step1_haplotyper/{}_sentieon.vcf".format(config['samples']['id']))
    
    if config['modules']['benchmarking']:
        run_all_files.append("Make_Vcf/step2_benchmarking/snp_compare/summary.txt")
        run_all_files.append("Make_Vcf/step2_benchmarking/indel_compare/summary.txt")

    if config['modules']['mate_pair_analysis']:
        run_all_files.append("Align/Distance_Between_Reads/mate_pairs_distance.pdf")

    if config['modules']['read_overlap_analysis']:
        run_all_files.append("Align/Distance_Between_Reads/distance_between_reads.pdf")
   
    if config['modules']['duplicate_plot']:
        run_all_files.append("Benchmarks/metrics.duplicate_plot.txt")

    if config['modules']['phasing']:
        for CHR in CHROMS:
            run_all_files.append("Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{}_hapblock_{}".format(config['samples']['id'], CHR))
            run_all_files.append("Make_Vcf/step4_longhap/LibOutDir/{}.vcf".format(CHR))

        run_all_files.extend(["Make_Vcf/step4_longhap/longhap_results.txt",
                              "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/hapcut_eval.txt",
                              "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/{}_hapcut.phased.vcf".format(config['samples']['id']),
                              "Make_Vcf/step4_longhap/{}_sentieon_pass_vars.vcf".format(config['samples']['id'])])

    return run_all_files


rule run_all:
    input:
        run_all_input
        

include: "/research/rv-02/home/eanderson/CGI_WGS_Pipeline/Snakefiles/calc_frag_len.snakefile"
include: "/research/rv-02/home/eanderson/CGI_WGS_Pipeline/Snakefiles/remove_sam_analysis.snakefile"
include: "/research/rv-02/home/eanderson/CGI_WGS_Pipeline/Snakefiles/make_vcf.snakefile"
include: "/research/rv-02/home/eanderson/CGI_WGS_Pipeline/Snakefiles/metrics.snakefile"
include: "/research/rv-02/home/eanderson/CGI_WGS_Pipeline/Snakefiles/stlfr.main.snakefile"
include: "/research/rv-02/home/eanderson/CGI_WGS_Pipeline/Snakefiles/splitreads.snakefile"
include: "/research/rv-02/home/eanderson/CGI_WGS_Pipeline/Snakefiles/phasing.snakefile"

shell.prefix('source /research/rv-02/home/eanderson/CGI_WGS_Pipeline/Data_and_Tools/bash_profile; ')


