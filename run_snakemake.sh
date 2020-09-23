#!/usr/bin/env bash

# Snakemake is installed in a virtual environment here
source /home/eanderson/Virtual_Envs/SnakeMake/bin/activate

# snakemake.err.txt usually has everything you need to figure out why something failed
# -j specifies the number of threads
# -k specifies that snakemake should continue to execute rules that are independant of any failures
# This saves time overall
# -s specifies the snakefile with the final targets
start=$(date +%s)
snakemake -j 110 -k -s /research/rv-02/home/eanderson/CGI_WGS_GATK_Pipeline/Snakefiles/run.all.snakefile 2> snakemake.err.txt

# Check that snakemake ran for at least an amount of time before ending
# This prevents a barrage of e-mails during testing or troubleshooting
if (((($(date +%s)-$start)/60) >= 20)); then
    echo "Snakemake has either finished or failed"
    # echo "Snakemake has either finished or failed" | \
    # Change the email address and uncomment
    # mailx -s "Snakemake Update" -a snakemake.err.txt "eaanderson@ucdavis.edu"
fi
