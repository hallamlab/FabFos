import json
import os, sys
from pathlib import Path
from Bio import SeqIO
from dataclasses import dataclass
import pandas as pd
from ..models import EndMappedContigs, RawContigs, EndSequences
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
    PIDENT_THRESH = 90
    for n, q in [
        ("forward", end_seqs.forward),
        ("reverse", end_seqs.reverse),
    ]:
        o = C.out_dir.joinpath(f"{n}.tsv")
        _cols = "qseqid sseqid nident qlen slen qstart qend sstart send".split(" ")
        os.system(f"""\
        echo "# {n}" >>{blast_log}          
        blastn \
            -perc_identity {PIDENT_THRESH} \
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
            if nident/qlen*100 < PIDENT_THRESH: continue
            front, back = hits.get(qseqid, (set(), set()))

            is_fwd = (sstart+send)/2 < slen/2
            collection = front if is_fwd else back
            collection.add(sseqid)
            
            hits[str(qseqid)] = front, back
        return hits
    
    ######################################
    # get contigs/scaffolds

    fhits = get_hits(hit_tables.forward)
    rhits = get_hits(hit_tables.reverse)
    with open(C.out_dir.joinpath("f.json"), "w") as j:
        json.dump({k: [list(a), list(b)] for k, (a, b) in fhits.items()}, j, indent=4)
    with open(C.out_dir.joinpath("r.json"), "w") as j:
        json.dump({k: [list(a), list(b)] for k, (a, b) in rhits.items()}, j, indent=4)
    
    end_ids = set(fhits)|set(rhits)

    stats = []
    mapped_file = Path("./mapped_contigs.fa")
    with open(mapped_file, "w") as f:
        def _wmapped(id, seq):
            f.write(f">{id}"+"\n")
            f.write(seq); f.write("\n")

        def _longest(contigs, flipped_contigs):
            def _it():
                for contig_id in contigs:
                    yield contig_id, all_contigs[contig_id]
                for contig_id in flipped_contigs:
                    yield contig_id, all_contigs[contig_id][::-1] # reverses string, flips contig back
            lid, longest, l = "", "", 0
            for cid, seq in _it():
                lseq = len(seq)
                if lseq > l:
                    l = lseq
                    longest = seq
                    lid = cid
            return lid, longest
        
        def _get_asm(contig_id):
            return  "_".join(contig_id.split("_")[1:])

        for i, id in enumerate(end_ids):
            fwd_front, fwd_back = fhits.get(id, (set(), set()))
            rev_front, rev_back = rhits.get(id, (set(), set()))
            matches = fwd_front&rev_back
            matches_contig_flipped = (rev_front&fwd_back) - matches # subtract matches just to be safe
            has_fwd_hits = len(fwd_front)+len(fwd_back) > 0
            has_rev_hits = len(rev_front)+len(rev_back) > 0

            if len(matches)+len(matches_contig_flipped) > 0:
                cid, longest = _longest(matches, matches_contig_flipped)
                quality, l, a = "full_match", len(longest), _get_asm(cid)
                _wmapped(f"{id} {quality} {a} length={l}", longest)
                
            elif has_fwd_hits and has_rev_hits:
                fid, fseq = _longest(fwd_front, fwd_back)
                rid, rseq = _longest(rev_back, rev_front)
                quality, l, a = "scaffolded", len(fseq)+len(rseq), f"{_get_asm(fid)},{_get_asm(rid)}"
                _wmapped(f"{id} {quality} {a} length>{l}", fseq+"-"+rseq)

            elif has_fwd_hits:
                cid, longest = _longest(fwd_front, fwd_back)
                quality, l, a = "forward_end_only", len(longest), _get_asm(cid)
                _wmapped(f"{id} {quality} {a} length>{l}", longest)

            elif has_rev_hits:
                cid, longest = _longest(rev_back, rev_front)
                quality, l, a = "reverse_end_only", len(longest), _get_asm(cid)
                _wmapped(f"{id} {quality} {a} length>{l}", longest)
            else:
                quality, l, a = "no_hits", None, None

            stats.append((id, quality, a, l))

        num_hits = len([i for i, q, l in stats if q != 'no_hits'])
        num_full_hits = len([i for i, q, l in stats if q == 'full_match'])
        C.log.info(f"found {num_hits} of {len(end_ids)}, {num_full_hits} were full matches")
        report_file = Path("./endmapping_report.csv")
        report = pd.DataFrame(stats, columns="id, mapping_quality, assembler, resolved_length".split(", "))
        report.to_csv(report_file, index=False)

        EndMappedContigs(mapped_file, report_file).Save(C.output)
