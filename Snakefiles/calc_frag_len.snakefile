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


def get_chroms(wildcards):
    if len(CHROMS) > 1:
        return ",".join(CHROMS)
    else:
        assert len(CHROMS) == 1
        return CHROMS


rule calc_frag_len:
    input:
        "Align/{}.sort.rmdup.bam".format(config['samples']['id'])
    output:
        "Calc_Frag_Length_{split}/frag_length_distribution.pdf",
        "Calc_Frag_Length_{split}/n_read_distribution.pdf",
        "Calc_Frag_Length_{split}/frag_and_bc_summary.txt",
        "Calc_Frag_Length_{split}/frags_per_bc.pdf"
    params:
        min_frag = config['calc_frag']['min_frag'],
        read_len = config['params']['read_len'],
        chroms = get_chroms,
        toolsdir = config['params']['toolsdir'],
        include_dups = config['calc_frag']['include_dups'],
        umi_analysis = config['modules']['umi_analysis']
    threads: 
        config['threads']['calc_frag']
    benchmark:
        "Benchmarks/calc_frag_len.calc_frag_len_{split}.txt"
    run:
        command = ["source /home/eanderson/.virtualenvs/General3/bin/activate ;",
                   "{params.toolsdir}/tools/calc_frag_len.py",
                   "--minfrag {params.min_frag}",
                   "--splitdist {wildcards.split}",
                   "--readlen {params.read_len}",
                   "--threads {threads}",
                   "--chroms {params.chroms}",
                   "--writeouttsvs ",
                   "--outdir Calc_Frag_Length_{wildcards.split}"]

        if params.include_dups:
            command.append("--includedups")

        command.append("{input}")
    
        shell(" ".join(command))
    
