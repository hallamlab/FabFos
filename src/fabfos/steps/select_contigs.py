import json
import os, sys
from pathlib import Path
from typing import Callable, Iterable
from Bio import SeqIO
from dataclasses import dataclass
import pandas as pd
from concurrent.futures import ProcessPoolExecutor as Exe
from ..models import ConsensusContigs, RawContigs, EndSequences
from .common import Init, FileSafeStr, Batchify

def _consensus(args):
    id, seqs, msa_temp, threads = args
    seqs_path = msa_temp.joinpath(f"{FileSafeStr(id)}.fa")
    msa_path = msa_temp.joinpath(f"{FileSafeStr(id)}.msa.clustal")
    cons_path = msa_temp.joinpath(f"{FileSafeStr(id)}.cons.fa")
    log_path = msa_temp.joinpath(f"{FileSafeStr(id)}.log")
    seqs = sorted(seqs, key=lambda t: len(t[1]), reverse=True)
    MAX_TO_COMPARE = 3
    with open(seqs_path, "w") as f:
        for i, (sid, s) in enumerate(seqs[:MAX_TO_COMPARE]):
            f.write(f">{sid}"+"\n")
            f.write(s); f.write("\n")
    os.system(f"""\
    mafft --thread {threads} --reorder --retree 2 --clustalout \
        {seqs_path} > {msa_path} 2>{log_path} \
    && cons -sequence {msa_path} -outseq {cons_path} >>{log_path} 2>&1
    """)
    return seqs_path

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
    db_folder = C.out_dir.joinpath("blast_dbs")
    os.makedirs(db_folder, exist_ok=True)
    blastdb = db_folder.joinpath("all_contigs_blastdb")
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
            # -evalue 1000 \
        os.system(f"""\
        echo "# {n}" >>{blast_log}          
        blastn \
            -perc_identity 90 \
            -num_threads {C.threads} \
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

    # _cols = "qseqid sseqid nident qlen slen".split(" ")
    # os.system(f"""\
    # echo "# pairwise" >>{blast_log}          
    # blastn \
    #     -evalue 1000 \
    #     -perc_identity 50 \
    #     -num_threads {C.threads} \
    #     -query {all_contigs_path} \
    #     -db {blastdb} \
    #     -outfmt "6 {' '.join(_cols)}" \
    #     -out {C.out_dir.joinpath("pairwise.tsv")} >>{blast_log} 2>&1
    # """)
    
    msa_temp = C.out_dir.joinpath("multiple_sequence_alignment")
    os.makedirs(msa_temp, exist_ok=True)

    fhits = get_hits(hit_tables.forward)
    rhits = get_hits(hit_tables.reverse)
    with open(C.out_dir.joinpath("f.json"), "w") as j:
        json.dump({k: [list(a), list(b)] for k, (a, b) in fhits.items()}, j, indent=4)
    with open(C.out_dir.joinpath("r.json"), "w") as j:
        json.dump({k: [list(a), list(b)] for k, (a, b) in rhits.items()}, j, indent=4)
    
    end_ids = set(fhits)|set(rhits)

    LOG_BUFFER = "      "
    THREADS_PER_JOB = 2
    n_parallel_jobs = C.threads//THREADS_PER_JOB

    def get_all_consensus():
        for i, id in enumerate(end_ids):
            if i > 3: break
            # print(f"{i+1} of {len(end_ids)}: {id} {LOG_BUFFER}", end="\r")
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
                yield id, _it()
            else:
                # partial_hits[id] = 
                # print(f"partial {id}      ")
                pass

    # if n_parallel_jobs>1:
    #     with Exe(max_workers=n_parallel_jobs) as exe:
    #         for res in exe.map(_consensus, iter((id, list(seqs), msa_temp, THREADS_PER_JOB) for id, seqs in get_all_consensus())):
    #             pass
    # else:
    #     for id, seqs in get_all_consensus():
    #         print(id)
    #         _consensus((id, seqs, msa_temp, C.threads))

    hits = {}
    partial_hits = {}



    # for k, (a, b) in hits.items():



    # return ConsensusContigs(*_hits.values())
    # contigs.Save(C.output)
