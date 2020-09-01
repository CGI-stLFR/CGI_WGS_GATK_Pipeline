rule make_hapcut_bam:
    input:
        "Align/{id}.sort.removedup_rm000.sam"
    output:
        "Align/{id}.sort.removedup_rm000.bam"
    shell:
        "samtools view -bh {input} > {output} && "
        "samtools index {output}"


rule split_hapcut_bam:
    input:
        "Align/{id}.sort.removedup_rm000.bam"
    output:
        "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.rmdup_{chr}.bam"
    params:
        id = config['samples']['id']
    shell:
        "samtools view -bh {input} {wildcards.chr} > "
            "{output} && "
            "samtools index {output}"


rule gzip_hapcut_vcf:
    input:
        "Make_Vcf/step1_haplotyper/{id}_sentieon_pass_vars.vcf"
    output:
        gz_vcf = "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars.vcf.gz",
        index = "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars.vcf.gz.tbi"
    params:
        bcftools = config['params']['bcftools']
    shell:
        "{params.bcftools} view -O z {input} > "
            "{output.gz_vcf} && "
            "tabix -p vcf -f {output.gz_vcf}"


rule split_hapcut_vcf:
    input:
        vcf = "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars.vcf.gz",
        index = "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars.vcf.gz.tbi"
    output:
        "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars_{chr}.vcf"
    params:
        bcftools = config['params']['bcftools'],
        id = config['samples']['id']
    shell:
        "{params.bcftools} view -O v "
            "-r {wildcards.chr} "
            "{input.vcf} > "
            "{output}"


rule get_hapcut_fragments:
    input:
        bam = "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.rmdup_{chr}.bam",
        vcf = "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars_{chr}.vcf"
    output:
        "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s1_unlinked_frag/{id}_unlinked_fragment_{chr}"
    params:
        toolsdir = config['params']['toolsdir']
    shell:
        "{params.toolsdir}/tools/HapCUT2/build/extractHAIRS --10X 1 "
            "--bam {input.bam} "
            "--VCF {input.vcf} "
            "--out {output}"


rule link_hapcut_fragments:
    input:
        bam = "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.rmdup_{chr}.bam",
        vcf = "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars_{chr}.vcf",
        frag = "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s1_unlinked_frag/{id}_unlinked_fragment_{chr}"
    output:
        "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s2_link_frag_files/{id}_linked_frag_{chr}"
    params:
        toolsdir = config['params']['toolsdir'],
        linkdist = config['params']['hapcut_link_dist']
    shell:
        "python3 {params.toolsdir}/tools/HapCUT2/utilities/LinkFragments.py --bam {input.bam} "
            "--VCF {input.vcf} "
            "--fragments {input.frag} "
            "--out {output} "
            "-d {params.linkdist}"


#def get_vcf(wildcards):
#    return "Make_Vcf/step3_hapcut/step2_split_vcf/{}_sentieon_pass_vars_{}.vcf".format(config['samples']['id'], wildcards.chr)


rule run_hapcut2:
    input:
        frag = "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s2_link_frag_files/{id}_linked_frag_{chr}",
        vcf = "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars_{chr}.vcf"
    output:
        blocks = "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{id}_hapblock_{chr}",
        vcf = "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{id}_hapblock_{chr}.phased.VCF"
    params:
        toolsdir = config['params']['toolsdir']
    shell:
        "{params.toolsdir}/tools/HapCUT2/build/HAPCUT2 --nf 1 "
            "--fragments {input.frag} "
            "--VCF {input.vcf} "
            "--output {output.blocks}"


rule evaluate_hapcut2:
    input:
        hapblock = expand("Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{id}_hapblock_{chr}", id=config['samples']['id'], chr=CHROMS),
        frag = expand("Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s2_link_frag_files/{id}_linked_frag_{chr}", id=config['samples']['id'], chr=CHROMS),
        vcf = expand("Make_Vcf/step3_hapcut/step2_split_vcf/{id}_sentieon_pass_vars_{chr}.vcf", id=config['samples']['id'], chr=CHROMS)
    output:
        "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/hapcut_eval.txt"
    params:
        toolsdir = config['params']['toolsdir'],
        truth_chroms = config['benchmark']['truth_vcf_dir'],
        ref_index = REF + '.fai'
    run:
        import sys
        from pathlib import Path

        if params.truth_chroms:
            truth_chroms = expand(params.truth_chroms + "/truth_chroms_{chr}.vcf", chr=CHROMS)
            shell("""python3 {toolsdir}/calculate_haplotype_statistics.py -h1 {input.hapblock} 
                     -v1 {input.vcf} -f1 {input.frag} 
                     -pv {truth_chroms} -c {params.ref_index} > {output}""")
        else:
            sys.path.append(str(Path(params.toolsdir) / 'tools'))
            from calculate_haplotype_statistics_no_truth import hapblock_vcf_error_rate_multiple

            err = hapblock_vcf_error_rate_multiple(input.hapblock, input.vcf, input.vcf, False)
            with open(output[0], "w") as outf:
                print(err, file=outf)


rule aggregate_hapcut_vcf:
    input:
        vcfs = expand("Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{id}_hapblock_{chr}.phased.VCF", id=config['samples']['id'], chr=CHROMS)
    output:
        agg_vcf = "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/{}_hapcut.phased.vcf".format(config['samples']['id'])
    run:
        import sys, re, glob
        from pathlib import Path

        contigs = []
        with open(input.vcfs[0], "r") as vcf_header, open(output.agg_vcf, "w") as outfile:
            for line in vcf_header:
                if line.startswith("#"):
                    print(line.strip(), file=outfile)
                    if "contig=" in line:
                        id_name = re.split("<|>", line.strip())[1]
                        keys = dict(parts.split("=") for parts in id_name.split(","))
                        contigs.append(keys['ID'])
                else:
                    break

        
        file_dict = {}
        for file in input.vcfs:
            name = Path(file).name
            tail = name.split("_hapblock_")[-1]
            contig = tail.split(".phased.VCF")[0]
            file_dict[contig] = file

        
        for contig in contigs:
            try:
                file = file_dict[contig]
            except KeyError:
                print(f"No file for contig {contig}, skipping...", file=sys.stderr)
                continue
            with open(file, "r") as vcf, open(output.agg_vcf, "a") as outfile:
                for line in vcf:
                    if line.startswith("#"):
                        continue
                    else:
                        print(line.strip(), file=outfile)


rule get_snps_for_longhap:
    input:
        "Make_Vcf/step1_haplotyper/{id}_sentieon_pass_vars.vcf"
    output:
        "Make_Vcf/step4_longhap/{id}_sentieon_pass_vars.vcf"
    params:
        bcftools = config['params']['bcftools']
    shell:
        "{params.bcftools} view -O v --type snps {input} > {output}"


rule generate_longhap_scripts:
    input:
        bam = "Align/{}.sort.rmdup.bam".format(config['samples']['id']),
        vcf = "Make_Vcf/step4_longhap/{}_sentieon_pass_vars.vcf".format(config['samples']['id'])
    output:
        "Make_Vcf/step4_longhap/LibShDir/run.{chr}.sh"
    params:
        win_size = config['longhap']['win_size'],
        max_len = config['longhap']['max_lfr_len'],
        min_bc = config['longhap']['min_barcode'],
        min_link = config['longhap']['min_link'],
        longhap = config['params']['longhap']
    run:
        from pathlib import Path
        import os

        root = Path("Make_Vcf/step4_longhap")
        dirs = ["LibShDir", "LibTmpDir", "LibOutDir"]

        for direc in dirs:
            direc_path = root/direc
            if not direc_path.exists():
                direc_path.mkdir()

        longhap_ex = Path(params.longhap)
        longhap_merge = longhap_ex.parent / "merge.pl"
        longhap_log2vcf = longhap_ex.parent / "log2vcf.pl"

        ref_index = REF + ".fai"
        ref_dict = {}
        with open(ref_index, "r") as ref_ind:
            for line in ref_ind:
                parts = line.strip().split()
                ref_dict[parts[0]] = int(parts[1])

        outfile = root / dirs[0] / "run.{}.sh".format(wildcards.chr)
        block_size = 10000000
        count_to = int(ref_dict[wildcards.chr])//block_size + 1
        with open(outfile, "w") as outf:

            for i in range(1, count_to):
                if i == 1:
                    beg = 1
                    end = block_size
                elif i < count_to:
                    beg = (i - 1) * block_size - 499999
                    end = i * block_size
                elif i == count_to:
                    beg = (i - 1) * block_size - 499999
                    end = int(ref_dict[wildcards.chr])

                print(f"perl {longhap_ex} {input.bam} {input.vcf} {wildcards.chr}:{beg}-{end} ",
                      f"{params.max_len} {params.min_bc} {params.min_link} ",
                      f"{root}/{dirs[1]}/{wildcards.chr}_{beg}_{end}.log 2> ",
                      f"{root}/{dirs[1]}/{wildcards.chr}_{beg}_{end}.error",
                      file=outf)

            print(f"perl {longhap_merge} {root}/{dirs[1]} {wildcards.chr} ",
                  f"{params.win_size} {root}/{dirs[2]}/{wildcards.chr}.log",
                  file=outf)
            print(f"cat {root/dirs[1]}/{wildcards.chr}_*.hete.barcodes | sort ",
                  f"| uniq > {root}/{dirs[2]}/{wildcards.chr}.hete.barcodes",
                  file=outf)
            print(f"cat {root/dirs[1]}/{wildcards.chr}_*.homo.barcodes | sort ",
                  f"| uniq > {root}/{dirs[2]}/{wildcards.chr}.homo.barcodes",
                  file=outf)
            print(f"cat {root/dirs[1]}/{wildcards.chr}_*.unknown.barcodes | sort",
                  f"| uniq >{root}/{dirs[2]}/{wildcards.chr}.unknown.barcodes",
                  file=outf)
            print(f"perl {longhap_log2vcf} {root}/{dirs[2]}/{wildcards.chr}.log ",
                  f"{root}/{dirs[2]}/{wildcards.chr}.hete.barcodes ",
                  f"{root}/{dirs[2]}/{wildcards.chr}.homo.barcodes ",
                  f"{input.vcf} > {root}/{dirs[2]}/{wildcards.chr}.vcf", file=outf)


rule run_longhap_shells:
    input:
        "Make_Vcf/step4_longhap/LibShDir/run.{chr}.sh"
    output:
        "Make_Vcf/step4_longhap/LibOutDir/{chr}.vcf"
    shell:
        "bash {input}"


rule evaluate_longhap:
    input:
        vcfs = expand("Make_Vcf/step4_longhap/LibOutDir/{chr}.vcf", chr=CHROMS),
        vcf = "Make_Vcf/step4_longhap/{}_sentieon_pass_vars.vcf".format(config['samples']['id'])
    output:
        "Make_Vcf/step4_longhap/longhap_results.txt"
    params:
        toolsdir = config['params']['toolsdir'],
        eval_vcfs = config['benchmark']['truth_vcf_dir'] # might be able to use the hapcut files.
    run:
        from pathlib import Path
        valid_pl = Path(params.toolsdir) / "tools" / "validate.pl"
        cal_ori = Path(params.toolsdir) / "tools" / "cal_ori_check.pl"
        eval_pl = Path(params.toolsdir) / "tools" / "evaluate.pl"

        commands = [['perl', str(valid_pl), input.vcf, 'Make_Vcf/step4_longhap/LibOutDir',
                     '|', 'tee', output[0]],
                    ['perl', str(cal_ori), input.vcf, 'Make_Vcf/step4_longhap/LibOutDir',
                     '|', 'tee', '-a', output[0]]]

        if params.eval_vcfs:
            commands.append(['perl', str(eval_pl), params.eval_vcfs,
                             'Make_Vcf/step4_longhap/LibOutDir', '|',
                             'tee', '-a', output[0]])

        for command in commands:
            shell(' '.join(command))

