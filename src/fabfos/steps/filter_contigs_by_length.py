import os, sys
from pathlib import Path
from Bio import SeqIO
import pandas as pd
from ..models import LenFilteredContigs, RawContigs
from .common import ClearTemp, Init

class Sequence:
    def __init__(self, id: str, seq: str, meta: dict|None=None) -> None:
        self.id = id
        self.seq = seq
        self.meta = meta

def Procedure(args):
    C = Init(args, __file__)
    raw_contigs = RawContigs.Load(C.NextArg())

    MIN_LEN = C.params.get("min_length", 1000)

    i = 1
    all_contigs: dict[str, Sequence] = {}
    all_contigs_path = C.out_dir.joinpath("all_contigs.fa")
    with open(all_contigs_path, "w") as f:
        for asm, con in raw_contigs.contigs.items():
            for e in SeqIO.parse(con, "fasta"):
                seq = str(e.seq)
                if len(seq)<MIN_LEN: continue
                k = f"{i:06}_{asm}"
                f.write(f">{k} length={len(seq)}"+"\n")
                f.write(seq); f.write("\n")
                all_contigs[k] = Sequence(k, seq, dict(assembler=asm))
                i += 1

    # make blastdb of contigs
    COLS = "qseqid, sseqid, nident, length, qlen, slen, qstart, qend, sstart, send".split(", ")
    db_folder = C.out_dir.joinpath("temp.blast_dbs")
    os.makedirs(db_folder, exist_ok=True)
    _blastdb = db_folder.joinpath("all_contigs")
    _duplicate_hits = C.out_dir.joinpath("pairwise.tsv")
    blast_log = C.out_dir.joinpath("blast_log.txt")
    C.shell(f"""\
    echo "# makeblastdb" >>{blast_log}
    makeblastdb \
        -dbtype nucl \
        -in {all_contigs_path} \
        -out {_blastdb} >>{blast_log} 2>&1
    echo "# blastn duplicaets" >>{blast_log}
    blastn \
        -perc_identity 50 \
        -num_threads {C.threads} \
        -query {all_contigs_path} \
        -db {_blastdb} \
        -outfmt "6 {' '.join(COLS)}" \
        -out {_duplicate_hits} >>{blast_log} 2>&1
    """)

    MAX_PIDENT = 99

    to_remove = set()
    try:
        similars = {}
        _df = pd.read_csv(_duplicate_hits, sep="\t", header=None)
        _df.columns = COLS
        for _, row in _df.iterrows():
            qseqid, sseqid, nident, length, qlen, slen, qstart, qend, sstart, send = row
            if qseqid == sseqid: continue
            if nident/qlen*100 >= MAX_PIDENT:
                similars[sseqid] = similars.get(sseqid, set())|{qseqid}

        for k, members in sorted(similars.items(), key=lambda t: len(t[1]), reverse=True):
            if k in to_remove: continue
            to_remove |= members

    except pd.errors.EmptyDataError:
        pass

    discard = open(C.out_dir.joinpath("discard.fna"), "w")
    contigs = open(C.root_workspace.joinpath("contigs.fna"), "w")
    for k, s in sorted(all_contigs.items(), key=lambda t: len(t[1].seq), reverse=True):
        f = discard if k in to_remove else contigs
        f.write(f">{k} length={len(s.seq)}"+"\n")
        f.write(s.seq); f.write("\n")
    discard.close()
    contigs.close()
    LenFilteredContigs(*[Path(fp.name) for fp in [contigs, discard]]).Save(C.expected_output)
    ClearTemp(C.out_dir)
