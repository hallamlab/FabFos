from __future__ import annotations
import os, sys
from typing import Any, Iterable
from pathlib import Path
from typing import Callable
from Bio import SeqIO
import pandas as pd
from ..models import EndMappedContigs, RawContigs, EndSequences
from .common import Init

class Sequence:
    def __init__(self, id: str, seq: str, meta: dict|None=None) -> None:
        self.id = id
        self.seq = seq
        self.meta = meta

class Hit:
    COLS = "qseqid, sseqid, nident, length, qlen, slen, qstart, qend, sstart, send".split(", ")

    def __init__(self, row: pd.Series, query_seqs: dict[str, Sequence], subject_seqs: dict[str, Sequence]) -> None:
        qseqid, sseqid, nident, length, qlen, slen, qstart, qend, sstart, send = row
        self._row = row
        self.qseqid = qseqid
        self.sseqid = sseqid
        self.length = length
        self.percent_query_identity = nident/qlen*100

        # to zero index
        qlen, slen, qstart, qend, sstart, send = [x-1 for x in [qlen, slen, qstart, qend, sstart, send]]
        
        self.query = query_seqs[qseqid]
        self.subject = subject_seqs[sseqid]
        # assert qstart < qend
        # if sstart<send:
        #     s, e = sstart-qstart, slen  # --->
        # else: # contig flipped
        #     s, e = 0, send+qstart       # <---
        # self.loc = s, e # in contig coordinates
        # self.is_left_to_right = not contig_flipped

    def _union(self):
        pass

    def Union(self, other_hits: Hit|Iterable[Hit]|None=None):
        if other_hits is None: other_hits = [] 
        elif isinstance(other_hits, Hit): other_hits = [other_hits]
        
    # check if hits are in opposite directions and facing each other
    def IsComplimentary(self, other: Hit):
        if self.is_left_to_right == (not other.is_left_to_right): return False # same direction
        s, e = self.Merge(other)
        if e < s: return False # not facing, <--- ---> instead of ---> -- <---
        return True
    
    def Merge(self, other: Hit|Iterable[Hit]):
        ss, se = self.loc
        if isinstance(other, Hit):
            os, oe = other.loc
        else:
            os = max(o.loc[0] for o in other)
            oe = min(o.loc[1] for o in other)
        return max(ss, os), min(se, oe)
        
def Procedure(args):
    C = Init(args)
    raw_contigs = RawContigs.Load(C.NextArg())
    end_seqs = EndSequences.Load(C.NextArg())
    assert end_seqs.forward is not None and end_seqs.reverse is not None

    ######################################
    # load end seqs and contigs as @Sequence

    MIN_LEN = 1000
    PIDENT_THRESH = 90

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

    all_ends: dict[str, tuple[Sequence, Sequence]] = {} # assumes all ends are paired
    _rev_ends = {}
    for e in SeqIO.parse(end_seqs.reverse, "fasta"):
        _rev_ends[e.id] = str(e.seq)
    for e in SeqIO.parse(end_seqs.forward, "fasta"):
        k = e.id
        all_ends[k] = Sequence(k, str(e.seq)), Sequence(k, _rev_ends[k])

    ######################################
    # blast ends against contigs

    # make blastdb of contigs
    db_folder = C.out_dir.joinpath("blast_dbs")
    os.makedirs(db_folder, exist_ok=True)
    _blastdb = db_folder.joinpath("all_contigs")
    blast_log = C.out_dir.joinpath("log.blast.txt")
    os.system(f"""\
    echo "# makeblastdb" >>{blast_log}
    makeblastdb \
        -dbtype nucl \
        -in {all_contigs_path} \
        -out {_blastdb} >>{blast_log} 2>&1
    """)

    # blast fwd and rev ends against contigs
    fwd_hits: dict[str, list[Hit]] = {}
    rev_hits: dict[str, list[Hit]] = {}
    for _name, _query_path, _hits in [
        ("hits_forward", end_seqs.forward, fwd_hits),
        ("hits_reverse", end_seqs.reverse, rev_hits),
    ]:
        o = C.out_dir.joinpath(f"{_name}.tsv")
        os.system(f"""\
        echo "# {_name}" >>{blast_log}
        blastn \
            -perc_identity {Hit.PIDENT_THRESH} \
            -num_threads {C.threads} \
            -query {_query_path} \
            -db {_blastdb} \
            -outfmt "6 {' '.join(Hit.COLS)}" \
            -out {o} >>{blast_log} 2>&1
        """)
        _df = pd.read_csv(o, sep="\t", header=None)
        _df.columns = Hit.COLS
        for _, row in _df.iterrows():
            h = Hit(row)
            if not h.valid: continue
            _hits[h.insert_id] = _hits.get(h.insert_id, [])+[h]
    
    ######################################
    # parse blast hits
    # defer scaffolding until later

    # helpers for next phase
    def _longest(contig_ids: Iterable[Sequence]):
        first = next(iter(contig_ids))
        longest, l = first, len(first.seq)
        for con in contig_ids:
            if len(con.seq) > l: longest, l = con, len(con.seq)
        return longest
    
    def _get_assembler(cid: str):
        meta = all_contigs[cid].meta
        UNK = "unk assembler"
        return meta.get("assembler", UNK) if meta is not None else UNK

    mapped_contigs = []
    _to_scaffold: list[tuple[str, dict[str, Hit], dict[str, Hit]]] = []
    for i, insert_id in enumerate(end_seqs.insert_ids):
        _dictify: Callable[[Any], dict[str, Hit]] = lambda hit_list: {h.contig_id:h for h in hit_list}
        fhits, rhits = _dictify(fwd_hits[insert_id]), _dictify(rev_hits[insert_id])

        # full contig match
        paired_hits = fhits.keys()&rhits.keys()
        if len(paired_hits)>0:
            def _get_contigs():
                for cid in paired_hits:
                    fwd, rev = fhits[cid], rhits[cid]
                    insert_id = fwd.insert_id
                    if not fwd.IsComplimentary(rev): continue
                    s, e = fwd.Merge(rev)
                    f, r = all_ends[insert_id]
                    s, e = s+len(f.seq), e-len(r.seq)
                    # merge ends and contig, chopping off where necessary
                    if s < e:
                        seq = f.seq+all_contigs[cid].seq[s:e]+r.seq[::-1]
                    else:
                        seq= f.seq+r.seq[::-1][e-s:]
                    yield Sequence(insert_id, seq, dict(
                        assembler=_get_assembler(cid),
                        quality = "full_match",
                    ))
            mapped_contigs.append(_longest(_get_contigs()))

        # fwd and rev hits, possible to scaffold
        elif len(fhits)>0 and len(rhits)>0:
            _to_scaffold.append((insert_id, fhits, rhits))

        # fwd hits only
        elif len(fhits)>0:
            def _get_fcontigs():
                for cid, hit in fhits.items():
                    s, e = hit.loc
                    end, _ = all_ends[hit.insert_id]
                    s += len(end.seq)
                    seq = end.seq+all_contigs[cid].seq[s:e]
                    yield Sequence(hit.insert_id, seq, dict(
                        assembler=_get_assembler(cid),
                        quality = "forward_end_only",
                    ))
            mapped_contigs.append(_longest(_get_fcontigs()))

        # rev hits only
        elif len(rhits)>0:
            def _get_rcontigs():
                for cid, hit in rhits.items():
                    s, e = hit.loc
                    _, end = all_ends[hit.insert_id]
                    e -= len(end.seq)
                    seq = all_contigs[cid].seq[s:e]+end.seq[::-1]
                    yield Sequence(hit.insert_id, seq, dict(
                        assembler=_get_assembler(cid),
                        quality = "reverse_end_only",
                    ))
            mapped_contigs.append(_longest(_get_rcontigs()))

    ######################################
    # blast fwd end hit contigs against
    # rev end hit contigs and scaffold

    SEP, MIN_SCAFFOLD_OVERLAP = "__", 200

    # make reverse end blastdb for scaffolding
    fwd_ends_for_scaffold = open(C.out_dir.joinpath("fwd_ends_for_scaffold.fa"), "w")
    rev_ends_for_scaffold = open(C.out_dir.joinpath("rev_ends_for_scaffold.fa"), "w")
    for insert_id, fhits, rhits in _to_scaffold:
        def _w(seq, cid, file):
            file.write(f">{insert_id}{SEP}{cid}"+"\n")
            file.write(seq);file.write("\n")
        fend, rend = all_ends[insert_id]
        fwd_list = [(fend, all_contigs[hit.contig_id], True, hit, fwd_ends_for_scaffold) for hit in fhits.values()]
        rev_list = [(rend, all_contigs[hit.contig_id], False, hit, rev_ends_for_scaffold) for hit in rhits.values()]
        for end, con, is_fwd, hit, file in  fwd_list+rev_list:
            s, e = hit.loc
            if is_fwd:
                seq = 
            _w(seq, should_flip, file)
    fwd_ends_for_scaffold.close(); rev_ends_for_scaffold.close()

    _blastdb = db_folder.joinpath("rev_end_for_scaffold")
    _scaffold_blast_results = C.out_dir.joinpath(f"scaffold_blast_results.tsv")
    _cols = "qseqid, sseqid, length, qlen, slen, qstart, qend, sstart, send".split(", ")
    os.system(f"""\
    echo "# scaffold makeblastdb" >>{blast_log}
    makeblastdb \
        -dbtype nucl \
        -in {rev_ends_for_scaffold.name} \
        -out {_blastdb} >>{blast_log} 2>&1
    echo "# scaffold" >>{blast_log}
    blastn \
        -perc_identity {Hit.PIDENT_THRESH} \
        -num_threads {C.threads} \
        -query {fwd_ends_for_scaffold.name} \
        -db {_blastdb} \
        -outfmt "6 {' '.join(_cols)}" \
        -out {_scaffold_blast_results} >>{blast_log} 2>&1
    """)
    _df = pd.read_csv(_scaffold_blast_results, sep="\t", header=None)
    _df.columns = _cols

    # get viable scaffolds if contigs overlap
    scaffolding_hits: dict[str, list[Sequence]] = {}
    for _, row in _df.iterrows():
        qseqid, sseqid, length, qlen, slen, qstart, qend, sstart, send = row
        if length < MIN_SCAFFOLD_OVERLAP: continue
        fwd_insert, fwd_contig, fflipped = qseqid.split(SEP)
        rev_insert, rev_contig, rflipped = sseqid.split(SEP)
        if fwd_insert != rev_insert: continue
        insert_id = fwd_insert
        fflipped, rflipped = fflipped=="True", rflipped=="True"
        fwd_l2r = (qstart<qend) != fflipped
        rev_l2r = (sstart<send) != rflipped
        if fwd_l2r != rev_l2r: continue # directionality based on ends should match

        def _get_seq(id, s, e, is_left):
            seq = all_contigs[id].seq
            if e<s:

            

        fwd_seq = all_contigs[fwd_contig].seq
        rev_seq = all_contigs[rev_contig].seq

        

    # ######################################
    # # get contigs/scaffolds

    # def _iterate_contigs(contigs, reverse_contigs):
    #     for contig_id in contigs:
    #         yield contig_id, all_contigs[contig_id]
    #     for contig_id in reverse_contigs:
    #         yield contig_id, all_contigs[contig_id][::-1] # reverses string, flips contig back
    # def _longest(contigs, reverse_contigs):
    #     lid, longest, l = "", "", 0
    #     for cid, seq in _iterate_contigs(contigs, reverse_contigs):
    #         if len(seq) > l: lid, longest, l = cid, seq, len(seq)
    #     return lid, longest
    
    # def _get_asm(contig_id):
    #     return  "_".join(contig_id.split("_")[1:])

    # mapped_file = open("./mapped_contigs.fa", "w")
    # def _wmapped(id, seq):
    #     mapped_file.write(f">{id}"+"\n")
    #     mapped_file.write(seq); mapped_file.write("\n")

    # fwd_scaf_file = open(C.out_dir.joinpath("scaffolding_fwd_ends.fa"), "w")
    # rev_scaf_file = open(C.out_dir.joinpath("scaffolding_rev_ends.fa"), "w")
    # scaffolds_to_check = []
    # def _wscaffold_check(id, fwd_front, fwd_back, rev_back, rev_front):
    #     def _w(cid, seq, f):
    #         f.write(f">{cid}"+"\n")
    #         f.write(seq); f.write("\n")
    #     for cid, seq in _iterate_contigs(fwd_front, fwd_back): _w(cid, seq, fwd_scaf_file)
    #     for cid, seq in _iterate_contigs(rev_back, rev_front): _w(cid, seq, rev_scaf_file)

    #     scaffolds_to_check.append((id, fwd_front, fwd_back, rev_back, rev_front))

    # stats = []
    # # parse mappings
    # for i, id in enumerate(end_ids):
    #     fwd_front, fwd_back = fhits.get(id, (set(), set()))
    #     rev_front, rev_back = rhits.get(id, (set(), set()))
    #     matches = fwd_front&rev_back
    #     matches_contig_flipped = (rev_front&fwd_back) - matches # subtract matches just to be safe
    #     has_fwd_hits = len(fwd_front)+len(fwd_back) > 0
    #     has_rev_hits = len(rev_front)+len(rev_back) > 0

    #     if len(matches)+len(matches_contig_flipped) > 0:
    #         cid, longest = _longest(matches, matches_contig_flipped)
    #         quality, l, a = "full_match", len(longest), _get_asm(cid)
    #         _wmapped(f"{id} {quality} {a} length={l}", longest)
            
    #     elif has_fwd_hits and has_rev_hits:
    #         # defer to later
    #         _wscaffold_check(id, fwd_front, fwd_back, rev_back, rev_front)
    #         continue

    #     elif has_fwd_hits:
    #         cid, longest = _longest(fwd_front, fwd_back)
    #         quality, l, a = "forward_end_only", len(longest), _get_asm(cid)
    #         _wmapped(f"{id} {quality} {a} length>{l}", longest)

    #     elif has_rev_hits:
    #         cid, longest = _longest(rev_back, rev_front)
    #         quality, l, a = "reverse_end_only", len(longest), _get_asm(cid)
    #         _wmapped(f"{id} {quality} {a} length>{l}", longest)
    #     else:
    #         quality, l, a = "no_hits", None, None

    #     stats.append((id, quality, a, l))

    # # find overlaps for scaffold or join with gap
    # _db = db_folder.joinpath("rev_end_hits")
    # _sca_hit_path = C.out_dir.joinpath("scaffolding_hits.tsv")
    # _cols = "qseqid sseqid nident qlen slen qstart qend sstart send".split(" ")
    # os.system(f"""\
    # echo "# scaffolding" >>{blast_log}
    # makeblastdb \
    #     -dbtype nucl \
    #     -in {rev_scaf_file.name} \
    #     -out {_db} >>{blast_log} 2>&1
    # blastn \
    #     -perc_identity 70 \
    #     -num_threads {C.threads} \
    #     -query {fwd_scaf_file.name} \
    #     -db {_db} \
    #     -outfmt "6 {' '.join(_cols)}" \
    #     -out {_sca_hit_path} >>{blast_log} 2>&1
    # """)
    # _df = pd.read_csv(_sca_hit_path, sep="\t", header=None)
    # _df.columns = _cols
    # _probable_rows = []
    # for i, row in _df.iterrows():
    #     qseqid, sseqid, nident, qlen, slen, qstart, qend, sstart, send = row
    #     # we expect the hit to be near the back of the query (fwd)
    #     # and near the front of the subject (rev)
    #     fwd_back = (qstart+qend)/2 > qlen/2
    #     rev_front = (sstart+send)/2 < slen/2
    #     if not fwd_back or not rev_front: continue
    #     # and that both hits are in the same direction
    #     fwd_flipped = qstart>qend
    #     rev_flipped = sstart>send
    #     if (not fwd_flipped and not rev_flipped) or (fwd_flipped and rev_flipped): continue
        
    #     _probable_rows.append(i)
    # _df = _df.iloc[_probable_rows]
    # _df.to_csv(_sca_hit_path, sep="\t", index=False)
        
    # for id, fwd_front, fwd_back, rev_back, rev_front in scaffolds_to_check:
    #     fwds = fwd_front|fwd_back
    #     revs = rev_front|rev_back

    # num_hits = len([i for i, q, a, l in stats if q != 'no_hits'])
    # num_full_hits = len([i for i, q, a, l in stats if q == 'full_match'])
    # C.log.info(f"found {num_hits} of {len(end_ids)}, {num_full_hits} were full matches")
    # report_file = Path("./endmapping_report.csv")
    # report = pd.DataFrame(stats, columns="id, mapping_quality, assembler, resolved_length".split(", "))
    # report.to_csv(report_file, index=False)

    # EndMappedContigs(mapped_file, report_file).Save(C.output)
