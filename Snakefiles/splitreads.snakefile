def get_fastqs_one(wildcards):
    from pathlib import Path
    fq_path = Path(config['samples']['fq_path'])
    lanes = config['samples']['lanes']
    fq_files = []
    for lane in lanes:
        loc_path = fq_path / lane
        for fp in loc_path.iterdir():
            if fp.is_file and str(fp).endswith("_1.fq.gz"):
                print(fp)
                fq_files.append(str(fp))
    return fq_files            


def get_fastqs_two(wildcards):
    from pathlib import Path
    fq_path = Path(config['samples']['fq_path'])
    lanes = config['samples']['lanes']
    fq_files = []
    for lane in lanes:
        loc_path = fq_path / lane
        for fp in loc_path.iterdir():
            if fp.is_file and str(fp).endswith("_2.fq.gz"):
                print(fp)
                fq_files.append(str(fp))
    return fq_files            


rule cat_read_one:
    input:
        get_fastqs_one
    output:
        "data/read_1.fq.gz"
    run:
        if len(config['samples']['lanes']) > 1:
            shell("cat {input} > {output}")
        else:
            shell("ln -s ../{input} {output}")


rule cat_read_two:
    input:
        get_fastqs_two
    output:
        "data/read_2.fq.gz"
    run:
        if len(config['samples']['lanes']) > 1:
            shell("cat {input} > {output}")
        else:
            shell("ln -s ../{input} {output}")


def determine_barcode_list(n_samples):
    import gzip
    import sys
    
    barcode_file = config['params']['toolsdir'] + config['params']['barcode']
    barcode_rc_file = config['params']['toolsdir'] + config['params']['barcode_RC']
    fastq_path = "data/read_2.fq.gz"
    
    def get_barcodes(barcodes_file):
        print(f"reading in {barcodes_file}", file=sys.stderr)
        nucs = ['A', 'C', 'G', 'T']
        barcodes = []
        with open(barcodes_file, "r") as bcs_file:
            for line in bcs_file:
                bc = line.strip().split()[0]
                for i in range(0, len(bc)):
                    for nuc in nucs:
                        bc_alt = list(bc)   
                        bc_alt[i] = nuc
                        barcodes.append("".join(bc_alt))

            barcodes_set = set(barcodes)
            return barcodes_set


    def get_fq_barcodes(fastq_path, n_samples):
        with gzip.open(fastq_path, "rb") as fastq:
            fq_bcs = []
            counter = 1
            print("Getting fastq barcodes", file=sys.stderr)
            for line in fastq:
                counter += 1
                if counter == n_samples + 1:
                    break
                if counter % 400000 == 0:
                    print(f"Read in {counter/4} lines", file=sys.stderr)
                if counter % 4 == 3:
                    seq = line.decode('utf-8').strip()
                    fq_bcs.append(seq[-10:])
                    fq_bcs.append(seq[-26:-16])
                    fq_bcs.append(seq[-42:-32])
                        
            return fq_bcs


    barcodes = get_barcodes(barcode_file)
    barcodes_rc = get_barcodes(barcode_rc_file)
    fq_bcs = get_fq_barcodes(fastq_path, n_samples)

    bcs_found = 0
    rc_bcs_found = 0
    for bc in fq_bcs:
        if bc in barcodes:
            bcs_found += 1
        if bc in barcodes_rc:
            rc_bcs_found += 1

    print(f"Barcodes: {bcs_found}\nRC_Barcodes: {rc_bcs_found}\n",file=sys.stderr)
    if bcs_found > rc_bcs_found:
        print("Barcode_list: {}".format(config['params']['barcode']), file=sys.stderr)
        return config['params']['barcode']
    elif rc_bcs_found > bcs_found:
        print("Barcode_list: {}".format(config['params']['barcode_RC']), file=sys.stderr)
        return config['params']['barcode_RC']
    else:
        if n_samples > 40000000:
            print("Can't determine correct barcode file for sample", file=sys.stderr)
            sys.exit(1)
        else:
            return determine_barcode_list(n_samples * 2)
            

def get_read_diff(read1, read2):
    import sys
    import gzip
    with gzip.open(read1, "r") as r1, gzip.open(read2, "r") as r2:
        r1.readline()
        r2.readline()
        seq1 = r1.readline().strip()
        seq2 = r2.readline().strip()
        if len(seq1) != config['params']['read_len']:
            print(f"interpreted read length ({len(seq1)}) does not match supplied read length ({config['params']['read_len']})", file=sys.stderr)
            sys.exit(1)
        if len(seq2) - len(seq1) == 42:
            print(f"42BP barcode detected", file=sys.stderr)
            return True
        elif len(seq2) - len(seq1) == 30:
            print(f"30BP barcode detected", file=sys.stderr)
            return False
        else:
            print(f"Unknown barcode length detected", file=sys.stderr)
            sys.exit(1)
        

rule split_reads:
    input:
        expand("data/read_{i}.fq.gz", i=range(1,3))
    output:
        expand("data/split_read.{i}.fq.gz", i=range(1,3))
    params:
        len = config['params']['read_len'],
        toolsdir = config['params']['toolsdir']
    run:
        params.barcode = determine_barcode_list(400000)
        if get_read_diff(input[0], input[1]):
            shell("perl {params.toolsdir}/tools/split_barcode_PEXXX_42_reads.pl "
                "{params.toolsdir}{params.barcode} "
                "{input} {params.len} data/split_read "
                "2> data/split_stat_read.err")
        else:
            shell("perl {params.toolsdir}/tools/split_barcode_PEXXX_30_reads.pl "
                "{params.toolsdir}{params.barcode} "
                "{input} {params.len} data/split_read "
                "2> data/split_stat_read.err")
