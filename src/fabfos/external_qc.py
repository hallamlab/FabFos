import os, sys
from pathlib import Path
import logging

def _read_counts(ws: Path, read_file: Path):
    #https://bioinformatics.stackexchange.com/questions/935/fast-way-to-count-number-of-reads-and-number-of-bases-in-a-fastq-file
    read_sizes = ws.joinpath("temp.readcount.txt")
    os.system(f"""\
        cat {read_file} \
        | awk 'NR % 4 == 2' \
        | wc -cl >{read_sizes} \
    """)

    num_reads, nucleotides = 0, 0
    with open(read_sizes) as f:
        toks = f.readline()[:-1].strip()
        if "\t" in toks: toks = toks.split("\t")
        else: toks = [t for t in toks.split(" ") if len(t)>0]
        num_reads, nucleotides = toks
    os.unlink(read_sizes)

    return int(nucleotides), int(num_reads)

# runs fastqc on reads
def Fastqc(ws: Path, reads: list[Path]):
    logging.info("Running Fastqc...")
    OUT = ws.joinpath("fastqc"); os.makedirs(OUT, exist_ok=True)
    LOG = OUT.joinpath("log.txt")
    os.system(f"""\
        fastqc -o {OUT} {" ".join(str(r) for r in reads)} 1>{LOG} 2>{LOG}
    """)
    logging.info("done.\n")
    return OUT

# gets read coverage of contigs and various assembly stats
def AssemblyStats(ws: Path, reads: list[Path], assembly: Path, cpus: int, paired_end=True):
    logging.info("Calculating assembly statistics...")
    LOG = "log.txt"
    OUT = ws.joinpath("temp_assembly_stats"); os.makedirs(OUT, exist_ok=True)
    READ_FILES = [r.absolute() for r in reads]
    ASM = assembly.absolute()
    nucleotides, read_count = 0, 0
    for r in READ_FILES:
        _n, _c = _read_counts(OUT, r)
        nucleotides+=_n
        read_count+=_c

    pe_params = "-x sr" if paired_end else ""

    COV = f"./assembly_coverage.tsv"
    COV_HEADER = "\t".join("contig, start, end, coverage".split(", "))
    os.system(f"""\
        cd {OUT}
        BAM=./temp.coverage.bam
        
        # align reads to assembly
        minimap2 -a {pe_params} --secondary=no {ASM} {' '.join([str(f) for f in READ_FILES])} 2>>{LOG} | samtools sort --threads {cpus} -o $BAM --write-index - 2>>{LOG}

        # read coverage
        echo "{COV_HEADER}">{COV}
        bedtools genomecov -ibam $BAM -bg >>{COV}
        mv {COV} ../
        samtools flagstat $BAM >./stats.txt

        # quast for asm stats
        quast -t {cpus} -o ./quast {ASM} 1>>{LOG}
        mv ./quast ../
    """)
    # with open(OUT.joinpath(COV)) as f:
    #     with open(ws.joinpath(COV), "w") as out:
    #         out.write("\t".join("contig, start, end, coverage".split(", ")))
    #         out.writelines(f.readlines())
    #         for l in f:
    #             # out.wr
    #             # con, s, e, cov = l.split("\t")
    #             # con = int(con.split("_")[-1])
    #             # con = f"{con:04}"
    #             # out.write("\t".join([con, s, e, cov]))
    # os.unlink(OUT.joinpath(COV))
        
    logging.info("done.\n")
    return OUT
