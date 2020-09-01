
CHROM_IDS = "Align/Distance_Between_Reads/" + REF.split('/').pop() + ".chrom.id"

rule remove_sam_duplicates:
    input:
        "Align/{id}.aln_mem.sam"
    output:
        "Align/{id}.aln.sam"
    params:
        readlen = config['params']['read_len']
    benchmark:
        "Benchmarks/remove_sam_analysis.remove_sam_duplicates.{id}.txt"
    run:
        shell("""
              perl -e '$readlen={params.readlen}; $a="";$b="";$id="";$n=0;
                while(<>){{if(!(/^@/)){{chomp;@t=split;
                if($t[0] ne $id){{if($n>=2){{print "$a\n$b\n";}}$a="$_"; $n=1;}}
                elsif(length($t[9]) == $readlen){{$b="$_"; $n++;}} $id=$t[0]; }}
                else{{print "$_";}}}}
                if($n>=2){{print "$a\n$b\n";}}' {input} > {output}
              """)


rule generate_chrom_ids:
    input:
        REF
    output:
        CHROM_IDS
    params:
        toolsdir = config['params']['toolsdir']
    benchmark:
        "Benchmarks/remove_sam_analysis.generate_chrom_ids.txt"
    shell:
        "perl {params.toolsdir}/tools/Fa_Stat.pl {input} {output}"


rule lfr_length_stats:
    input:
        chroms = CHROM_IDS,
        sam = "Align/{}.aln.sam".format(config['samples']['id'])
    output:
        "Align/Distance_Between_Reads/reads_on_ref.pos"
    params:
        toolsdir = config['params']['toolsdir'],
        readlen = config['params']['read_len']
    benchmark:
        "Benchmarks/remove_sam_analysis.lfr_length_stats.txt"
    shell:
        "perl {params.toolsdir}/tools/LFR_Length_Stat_new5.pl "
            "{input.chroms} "
            "{input.sam} "
            "0 {params.readlen} none "
            "{output}"


rule overlap_reads_plot:
    input:
        "Align/Distance_Between_Reads/reads_on_ref.pos"
    output:
        "Align/Distance_Between_Reads/overlap_reads.plot"
    params:
        toolsdir = config['params']['toolsdir']
    benchmark:
        "Benchmarks/remove_sam_analysis.overlap_reads_plot.txt"
    shell:
        "perl {params.toolsdir}/tools/overlap_reads.pl {input} {output}"


rule overlap_hist:
    input:
        "Align/Distance_Between_Reads/overlap_reads.plot"
    output:
        "Align/Distance_Between_Reads/overlap_reads.plot.R"
    params:
        toolsdir = config['params']['toolsdir']
    benchmark:
        "Benchmarks/remove_sam_analysis.overlap_hist.txt"
    shell:
        "perl {params.toolsdir}/tools/Histogram_plot.pl "
            "{input} 0 1 1000 {output}"


rule reads_on_ref_count:
    input:
        "Align/Distance_Between_Reads/reads_on_ref.pos"
    output:
        "Align/Distance_Between_Reads/reads_on_ref.pos.num"
    benchmark:
        "Benchmarks/remove_sam_analysis.reads_on_ref_count.txt"
    shell:
        "wc -l {input} > {output}"


rule overlap_reads_plot_2:
    input:
        count = "Align/Distance_Between_Reads/reads_on_ref.pos.num",
        plot = "Align/Distance_Between_Reads/overlap_reads.plot.R"
    output:
        "Align/Distance_Between_Reads/overlap_reads.plot2.R"
    params:
        readlen = config['params']['read_len']
    benchmark:
        "Benchmarks/remove_sam_analysis.overlap_reads_plot_2.txt"
    run:
        shell("""
              perl -e '$readlen={params.readlen}; open IN1,"{input.count}";
                $N=<IN1>;@t=split(/\s+/,$N);$num=$t[0]; close(IN1);
                open IN2,"{input.plot}";
                while(<IN2>){{chomp;@t=split;
                if(/^0/){{$s0=$t[1];}}
                else{{$x=$t[0]-$readlen;
                $r=100*$t[1]/($num-$s0);
                print "$x\t$t[1]\t$r\n"; }} }} close(IN2);' > {output}
              """)


rule plot_read_overlap:
    input:
        "Align/Distance_Between_Reads/overlap_reads.plot2.R"
    output:
        "Align/Distance_Between_Reads/distance_between_reads.pdf"
    params:
        toolsdir = config['params']['toolsdir']
    benchmark:
        "Benchmarks/remove_sam_analysis.plot_read_overlap.txt"
    shell:
        "{params.toolsdir}/tools/R-3.2.1/bin/Rscript "
            "{params.toolsdir}/tools/Distance_between_reads.R "
            "{input} 8 {output}"


rule mate_pair_distance:
    input:
        chroms = CHROM_IDS,
        sam = "Align/{}.aln.sam".format(config['samples']['id'])
    output:
        "Align/Distance_Between_Reads/mate_pairs_distance.plot"
    params:
        toolsdir = config['params']['toolsdir'],
        readlen = config['params']['read_len']
    benchmark:
        "Benchmarks/remove_sam_analysis.mate_pair_distance.txt"
    shell:
        "perl {params.toolsdir}/tools/mate_pairs_distance.pl "
            "{input.chroms} "
            "{input.sam} "
            "{params.readlen} "
            "Align/Distance_Between_Reads/mate_pairs_distance"


rule mate_pair_hist:
    input:
        "Align/Distance_Between_Reads/mate_pairs_distance.plot"
    output:
        "Align/Distance_Between_Reads/mate_pairs_distance.plot.R"
    params:
        toolsdir = config['params']['toolsdir']
    benchmark:
        "Benchmarks/remove_sam_analysis.mate_pair_hist.txt"
    shell:
        "perl {params.toolsdir}/tools/Histogram_plot.pl "
            "{input} 0 1000 300000000 {output}"


rule mate_pair_count:
    input:
        "Align/Distance_Between_Reads/mate_pairs_distance.plot"
    output:
        "Align/Distance_Between_Reads/mate_pairs_distance.num"
    benchmark:
        "Benchmarks/remove_sam_analysis.mate_pair_count.txt"
    shell:
        "wc -l {input} > {output}"


rule mate_pair_plot_2:
    input:
        count = "Align/Distance_Between_Reads/mate_pairs_distance.num",
        plot = "Align/Distance_Between_Reads/mate_pairs_distance.plot.R"
    output:
        "Align/Distance_Between_Reads/mate_pairs_distance.plot2.R"
    params:
        readlen = config['params']['read_len']
    benchmark:
        "Benchmarks/remove_sam_analysis.mate_pair_plot_2.txt"
    run:
        shell("""
              perl -e '$readlen={params.readlen}; open IN1,"{input.count}";
                $N=<IN1>;@t=split(/\s+/,$N);$num=$t[0]; close(IN1);
                open IN2,"{input.plot}";
                while(<IN2>){{chomp;@t=split;
                if(/^0/){{$s0=$t[1];}}
                else{{$x=$t[0]-$readlen;
                $r=100*$t[1]/($num-$s0);
                print "$x\t$t[1]\t$r\n"; }} }} close(IN2);' > {output}
              """)


rule plot_mate_pairs:
    input:
        "Align/Distance_Between_Reads/mate_pairs_distance.plot2.R"
    output:
        "Align/Distance_Between_Reads/mate_pairs_distance.pdf"
    params:
        toolsdir = config['params']['toolsdir']
    benchmark:
        "Benchmarks/remove_sam_analysis.plot_mate_pairs.txt"
    shell:
        "{params.toolsdir}/tools/R-3.2.1/bin/Rscript "
            "{params.toolsdir}/tools/Mate_pairs_distance.R  "
            "{input} 300000000 10 {output}"
