import json
import os, sys
from pathlib import Path
from typing import Iterable
from Bio import SeqIO
from dataclasses import dataclass
import pandas as pd
from ..models import ConsensusContigs, RawContigs, EndSequences
from .common import Init

def Procedure(args):
    C = Init(args)
    raw_contigs = RawContigs.Load(C.NextArg())
    end_seqs = EndSequences.Load(C.NextArg())

    i = 1
    all_contigs = {}
    all_contigs_path = C.out_dir.joinpath("all_contigs.fa")
    with open(all_contigs_path, "w") as f:
        def _w(k, seq):
            seq = str(seq)
            f.write(f">{k} length={len(seq)}"+"\n")
            f.write(seq)
            f.write("\n")
            all_contigs[k] = seq

        for asm, con in raw_contigs.contigs.items():
            for e in SeqIO.parse(con, "fasta"):
                _w(f"{i:06}_{asm}", e.seq)
                i += 1

    assert end_seqs.forward is not None and end_seqs.reverse is not None

    ######################################
    # map end seqs
    blastdb = C.out_dir.joinpath("all_contigs_blastdb")
    blast_log = C.out_dir.joinpath("log.blast.txt")
    os.system(f"""\
    echo "# makeblastdb" >>{blast_log}
    makeblastdb \
        -dbtype nucl \
        -in {all_contigs_path} \
        -out {blastdb} >>{blast_log} 2>&1
    """)

    @dataclass
    class Hits:
        forward: pd.DataFrame
        reverse: pd.DataFrame
    _hits = {}
    for n, q in [
        ("forward", end_seqs.forward),
        ("reverse", end_seqs.reverse),
    ]:
        o = C.out_dir.joinpath(f"{n}.tsv")
        _cols = "qseqid sseqid nident qlen slen qstart qend sstart send".split(" ")
        os.system(f"""\
        echo "# {n}" >>{blast_log}          
        blastn \
            -evalue 1000 \
            -perc_identity 50 \
            -query {q} \
            -db {blastdb} \
            -outfmt "6 {' '.join(_cols)}" \
            -out {o} >>{blast_log} 2>&1
        """)
        _df = pd.read_csv(o, sep="\t", header=None)
        _df.columns = _cols
        _hits[n] = _df
    hit_tables = Hits(**_hits)

    def get_hits(df: pd.DataFrame):
        hits = {}
        for _, row in df.iterrows():
            qseqid, sseqid, nident, qlen, slen, qstart, qend, sstart, send = row
            front, back = hits.get(qseqid, (set(), set()))

            is_fwd = (sstart+send)/2 < slen/2
            collection = front if is_fwd else back
            collection.add(sseqid)
            
            hits[str(qseqid)] = front, back
        return hits
    
    # hits = get_hits(raw_hits.forward, is_fwd=True)
    # for k, (f, r) in get_hits(raw_hits.reverse, is_fwd=False).items():
    #     prevf, prevr = hits.get(k, (set(), set()))
    #     hits[k] = prevf|f, prevr|r
    
    ######################################
    # MSA and consensus contigs

    _cols = "qseqid sseqid nident qlen slen".split(" ")
    os.system(f"""\
    echo "# pairwise" >>{blast_log}          
    blastn \
        -evalue 1000 \
        -perc_identity 50 \
        -query {all_contigs_path} \
        -db {blastdb} \
        -outfmt "6 {' '.join(_cols)}" \
        -out {C.out_dir.joinpath("pairwise.tsv")} >>{blast_log} 2>&1
    """)
    
    msa_temp = C.out_dir.joinpath("multiple_sequence_alignment")
    os.makedirs(msa_temp, exist_ok=True)
    def file_safe(s: str):
        WL = "-_."
        return "".join(ch if (ch.isalnum() or ch in WL) else "_" for ch in s)
    def consensus(id, seqs: Iterable):
        seqs_path = msa_temp.joinpath(f"{file_safe(id)}.fa")
        with open(seqs_path, "w") as f:
            for i, (sid, s) in enumerate(seqs):
                # f.write(f">{i:03}"+"\n")
                f.write(f">{sid}"+"\n")
                f.write(s); f.write("\n")
        # echo "# pairwise" >>{blast_log}
        _cols = "qseqid sseqid nident qlen slen".split(" ")
        _db = msa_temp.joinpath(f"{file_safe(id)}_db")
        _out = msa_temp.joinpath(f"{file_safe(id)}_aln.tsv")
        os.system(f"""\
        makeblastdb \
            -dbtype nucl \
            -in {seqs_path} \
            -out {_db} >/dev/null 2>&1
        blastn \
            -evalue 1000 \
            -perc_identity 50 \
            -query {all_contigs_path} \
            -db {_db} \
            -outfmt "6 {' '.join(_cols)}" \
            -out {_out}
        """)
        # cons_path = msa_temp.joinpath(f"{file_safe(id)}.consensus.fa")
        # os.system(f"""\
        # mafft --thread {C.threads} --retree 1 --maxiterate 2 \
        #     {seqs_path} > {cons_path}
        # """)

    fhits = get_hits(hit_tables.forward)
    rhits = get_hits(hit_tables.reverse)
    with open(C.out_dir.joinpath("f.json"), "w") as j:
        json.dump({k: [list(a), list(b)] for k, (a, b) in fhits.items()}, j, indent=4)
    with open(C.out_dir.joinpath("r.json"), "w") as j:
        json.dump({k: [list(a), list(b)] for k, (a, b) in rhits.items()}, j, indent=4)
    
    hits = {}
    partial_hits = {}
    end_ids = set(fhits)|set(rhits)
    for id in end_ids:
        fwd_front, fwd_back = fhits.get(id, (set(), set()))
        rev_front, rev_back = rhits.get(id, (set(), set()))
        fwd_matches = fwd_front&rev_back
        rev_matches = (rev_front&fwd_back) - fwd_matches # just to be safe

        if len(fwd_matches)+len(rev_matches) > 0:
            def _it():
                for contig_id in fwd_matches:
                    yield contig_id, all_contigs[contig_id]
                for contig_id in rev_matches:
                    yield contig_id, all_contigs[contig_id][::-1] # reverses string
            consensus(id, _it())
        else:
            # partial_hits[id] = 
            pass


    # for k, (a, b) in hits.items():



    # return ConsensusContigs(*_hits.values())
    # contigs.Save(C.output)
