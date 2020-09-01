#!/usr/bin/env bash

source /home/eanderson/Virtual_Envs/SnakeMake/bin/activate

start=$(date +%s)
snakemake -j 110 -k -s /research/rv-02/home/eanderson/CGI_WGS_Pipeline/Snakefiles/run.all.snakefile 2> snakemake.err.txt


if (((($(date +%s)-$start)/60) >= 20)); then
    echo "Snakemake has either finished or failed" | \
    mailx -s "Snakemake Update" -a snakemake.err.txt "eaanderson@ucdavis.edu"
fi
