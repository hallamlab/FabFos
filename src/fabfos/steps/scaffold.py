from __future__ import annotations
import os, sys
from dataclasses import dataclass
from typing import Any, Iterable, overload
from enum import Enum, auto
from pathlib import Path
from typing import Callable
from Bio import SeqIO
from Bio.Seq import Seq as RawSeq
import numpy as np
import pandas as pd
from ..models import Scaffolds, RawContigs, EndSequences
from .common import ClearTemp, Init

class Sequence:
    def __init__(self, id: str, seq: RawSeq, meta: dict|None=None) -> None:
        self.id = id
        self.raw_seq = seq
        self.meta = meta

    def __len__(self):
        return len(self.raw_seq)

    def Forward(self) -> RawSeq:
        return self.raw_seq
    
    def ReverseCompliment(self) -> RawSeq:
        return self.raw_seq.reverse_complement(inplace=False)

class Hit:
    COLS = "qseqid, sseqid, nident, length, qlen, slen, qstart, qend, sstart, send".split(", ")

    def __init__(self, row: pd.Series, query_seqs: dict[str, Sequence], subject_seqs: dict[str, Sequence]) -> None:
        qseqid, sseqid, nident, length, qlen, slen, qstart, qend, sstart, send = row
        self._row = row
        self.qseqid = qseqid
        self.sseqid = sseqid
        self.length = length
        self.percent_query_mapped = nident/qlen*100

        # convert to zero index open interval coordinates
        qstart-=1
        if sstart < send:
            sstart-=1
        else:
            send-=1
        
        # precalculate truncation of subject by query
        if qstart > qend or sstart > send: # one flipped
            s, e = sstart+qstart, 0 # sstart > send and so is further right
        else:
            s, e = sstart-qstart, slen
        self._subject_trunc_loc = s, e

        # precalculate union, uses query seq for overlap
        if qstart > qend or sstart > send:
            self._union_loc = 0, qend, sstart, 0 # qa1, qb1, sa1, sb1
        else:
            self._union_loc = None # must be flipped, since contigs oriented with ends at start

        self.query = query_seqs[qseqid]
        self.subject = subject_seqs[sseqid]

    def _cut(self, seq: Sequence, start: int, end: int):
        if start > end:
            return seq.ReverseCompliment()[len(seq)-start:len(seq)-end]
        else:
            return seq.Forward()[start:end]

    def _get_edges_and_overflows(self, s:int, e:int):
        if s < e:
            left_overflow = abs(min(s, 0))
            right_overflow = max(len(self.subject), e) - len(self.subject)
            s += left_overflow
            e -= right_overflow
        else:
            left_overflow = max(len(self.subject), s) - len(self.subject)
            right_overflow = abs(min(e, 0))
            s -= left_overflow
            e += right_overflow

        return s, e, left_overflow, right_overflow

    # when both hits mapped to the same contig
    def SeqFromPairedHit(self, other: Hit) -> RawSeq:
        # orienting self as forward end hit
        fs, fe = self._subject_trunc_loc
        rs, re = other._subject_trunc_loc
        # query should be the end seqs and subject is never flipped
        # |(re)<--------(fs)-->>(fwd end)----------(rev end)<--(rs)-------->>(fe)|
        # self.subject should == other.subject
        # so final scaffold is: (fwd end) + (assembly truncated to where the ends mapped) + (rev end)

        if (fs > fe) == (rs > re): # both hit are in the same direction (== is !XOR)
            return RawSeq("") # end seqs should be facing; this is a false positive
            # return self.query.Forward()+self._cut(self.subject, fs, rs)+other.query.Forward() # return anyway for testing
        if fs < fe and fs > rs: # fwd is left but rev is further left
            return RawSeq("")
        if rs < re and rs > fs: # rev is left but fwd is further left
            return RawSeq("")
        
        s, e, left_overflow, right_overflow = self._get_edges_and_overflows(fs, rs) # _cut will handle reverse compliment if rs < fs
        if left_overflow>0:
            left_pad = self._cut(self.query, 0, left_overflow)
        else:
            left_pad = RawSeq("")
        if right_overflow>0:
            other_q = other.query
            right_pad_rc = self._cut(other_q, right_overflow, 0) # reverse compliments
        else:
            right_pad_rc = RawSeq("")

        # if self.qseqid=="case_facing_but_past":
        #     print(self.qseqid, s, e, "|", fs, fe, rs, re)
        #     print(self.qseqid, left_overflow, right_overflow)
        #     print(self.qseqid, len(left_pad), len(self._cut(self.subject, s, e)), len(right_pad_rc))
        return left_pad+self._cut(self.subject, s, e)+right_pad_rc

    def AlignSubjectToQuery(self) -> RawSeq:
        s, e, left_overflow, right_overflow = self._get_edges_and_overflows(*self._subject_trunc_loc)
        if left_overflow>0:
            left_pad = self._cut(self.query, 0, left_overflow)
        else:
            left_pad = RawSeq("")

        if right_overflow>0:
            right_pad = self._cut(self.query, len(self.query)-right_overflow, len(self.query))
        else:
            right_pad = RawSeq("")
        if self.qseqid=="case_fwd_only": 
            print(self.qseqid, s, e, "|", left_overflow, right_overflow, len(self.subject))
            print(self.qseqid, len(left_pad), len(self._cut(self.subject, s, e)), len(right_pad))
            # print(self.qseqid, self.subject.Forward()[s:e])
            # print(self.qseqid, self._cut(self.subject, s, e))
        return left_pad+self._cut(self.subject, s, e)+right_pad

    def Union(self):
        if self._union_loc is None: return RawSeq("")
        qa, qb, sa, sb = self._union_loc
        overlap = self.length
        qseq = self._cut(self.query, qa, qb)
        sseq = self._cut(self.subject, sa, sb) # sa > sb, reverse compliment
        # print(self.qseqid, qa, qb, sa, sb, self.length)
        # print(self.qseqid, len(qseq), len(sseq))
        return qseq+sseq[overlap:]
        
def Procedure(args):
    C = Init(args, __file__)
    raw_contigs = RawContigs.Load(C.NextArg())
    end_seqs = EndSequences.Load(C.NextArg())
    expected_insert_length = C.params["exp_length"]
    expected_insert_length_range = C.params["exp_length_range"]
    min_contig_length = C.params["min_contig_length"]
    GAP_CHAR = C.params["gap_str"]
    ends_facing = C.params.get("ends_facing", False)
    PIDENT_THRESH = C.params.get("pident", 90)
    assert end_seqs.forward is not None and end_seqs.reverse is not None

    ######################################
    # ensure ends are facing

    if not ends_facing:
        new_rev_path = C.out_dir.joinpath("reversed_end_seqs.reverse_compliment.fna")
        with open(new_rev_path, "w") as _out:
            for seq in SeqIO.parse(end_seqs.reverse, "fasta"):
                _out.write(f">{seq.description}\n")
                _out.write(str(seq.seq.reverse_complement(inplace=False))+"\n")
        end_seqs.reverse = new_rev_path

    ######################################
    # load end seqs and contigs as @Sequence

    all_contigs: dict[str, Sequence] = {}
    original_contig_id = {}
    all_contigs_path = C.out_dir.joinpath("all_contigs.fa")
    temp_counts_dir = C.out_dir/"temp_counts"
    if temp_counts_dir.exists(): C.shell(f"rm -rf {temp_counts_dir}")
    os.makedirs(temp_counts_dir, exist_ok=True)
    i = 1
    with open(all_contigs_path, "w") as f:
        count = 0 # to get number of digits for id of contigs
        for asm, file_path in raw_contigs.contigs.items():
            _temp = temp_counts_dir/f"{asm}.count"
            C.shell(f"cat {file_path} | grep \">\" | wc -l > {_temp}")
            with open(_temp, "r") as _count_file:
                count += int(_count_file.read())
        digits = int(np.ceil(np.log10(count)))

        for asm, file_path in raw_contigs.contigs.items():
            for e in SeqIO.parse(file_path, "fasta"):
                if len(e.seq)<min_contig_length: continue
                k = f"C{i:0{digits}}"
                original_contig_id[k] = e.id
                f.write(f">{k} asm={asm} length={len(e.seq)} | {e.description}"+"\n")
                f.write(str(e.seq)); f.write("\n")
                all_contigs[k] = Sequence(k, e.seq, dict(assembler=asm))
                i += 1

    fwd_ends: dict[str, Sequence] = {}
    rev_ends: dict[str, Sequence] = {}
    for f, d in [
        (end_seqs.forward, fwd_ends),
        (end_seqs.reverse, rev_ends),
    ]:
        for e in SeqIO.parse(f, "fasta"):
            k = str(e.id)
            d[k] = Sequence(k, e.seq)

    ######################################
    # blast ends against contigs

    db_folder = C.out_dir.joinpath("temp.blast_dbs")
    os.makedirs(db_folder, exist_ok=True)
    blast_log = C.out_dir.joinpath("blast_log.txt")
    C.log.info(f"blasting [{len(end_seqs.insert_ids)}] ends against [{len(all_contigs)}] contigs")

    # make blastdb of contigs
    _blastdb = db_folder.joinpath("all_contigs")
    C.shell(f"""\
    echo "# makeblastdb" >>{blast_log}
    makeblastdb \
        -dbtype nucl \
        -in {all_contigs_path} \
        -out {_blastdb} >>{blast_log} 2>&1
    """)

    # blast fwd and rev ends against contigs
    num_hits = 0
    fwd_hits: dict[str, list[Hit]] = {}
    rev_hits: dict[str, list[Hit]] = {}
    for _name, _query_path, _hits, _ends in [
        ("hits_forward", end_seqs.forward, fwd_hits, fwd_ends),
        ("hits_reverse", end_seqs.reverse, rev_hits, rev_ends),
    ]:
        o = C.out_dir.joinpath(f"{_name}.tsv")
        C.shell(f"""\
        echo "# {_name}" >>{blast_log}
        blastn \
            -evalue 999 \
            -perc_identity {PIDENT_THRESH} \
            -num_threads {C.threads} \
            -query {_query_path} \
            -db {_blastdb} \
            -outfmt "6 {' '.join(Hit.COLS)}" \
            -out {o} >>{blast_log} 2>&1
        """)
        try:
            _df = pd.read_csv(o, sep="\t", header=None)
            _df.columns = Hit.COLS
            for _, row in _df.iterrows():
                h = Hit(row, _ends, all_contigs)
                if h.percent_query_mapped < PIDENT_THRESH:
                    # print(h.qseqid, h.percent_query_mapped)
                    continue
                k = h.qseqid
                _hits[k] = _hits.get(k, [])+[h]
                num_hits += 1
        except pd.errors.EmptyDataError:
            pass
        
    ######################################
    # parse blast hits
    # defer scaffolding until later

    # helpers for next phase
    def _longest(sequences: Iterable[Sequence]):
        first = next(iter(sequences))
        def _get(x: Sequence):
            return x.Forward()
        longest, l = first, len(_get(first))
        for con in sequences:
            if len(_get(con)) > l: longest, l = con, len(_get(con))
        return longest
    
    def _get_assembler(cid: str|None):
        UNK = "unk_assembler"
        if cid is None: return UNK
        meta = all_contigs[cid].meta
        return meta.get("assembler", UNK) if meta is not None else UNK

    # score multipliers for choosing "longest" resolved contigs / scaffolds
    class QUALITY(float, Enum):
        FULL_MATCH =        10
        FULL_SCAFFOLD =     7
        GAPPED_SCAFFOLD =   5
        FORWARD_END_ONLY =  3.001 # enums can't have same value...
        REVERSE_END_ONLY =  3
        NO_HITS =           0

    C.log.info(f"mapping [{num_hits}] hits ")
    mapped_contigs: dict[str, list[Sequence]] = {}
    _to_scaffold: list[tuple[str, list[tuple[str, Hit]], list[tuple[str, Hit]]]] = []
    for i, insert_id in enumerate(end_seqs.insert_ids):
        def _dictify(hit_dict: dict[str, list[Hit]]):
            _hits: dict[str, list[Hit]] = {}
            for h in hit_dict.get(insert_id, []):
                _hits[h.sseqid] = _hits.get(h.sseqid, [])+[h]
            return _hits
        fhits, rhits = _dictify(fwd_hits), _dictify(rev_hits)

        # if insert_id == "CEC_1723":
        #     print([(h.sseqid, h.) for h in fhits.values()])

        # full contig match
        _candidates = []
        paired_hits = fhits.keys()&rhits.keys()
        if len(paired_hits)>0:
            def _iterate_hit_pairs(cid: str):
                for fwd in fhits[cid]:
                    for rev in rhits[cid]:
                        yield fwd, rev
            def _get_contigs():
                for cid in paired_hits:
                    for fwd, rev in _iterate_hit_pairs(cid):
                        seq = fwd.SeqFromPairedHit(rev)
                        # if insert_id == "CEC_1723":
                        #     print(len(seq), cid, fwd._subject_trunc_loc, rev._subject_trunc_loc)
                        yield Sequence(insert_id, seq, dict(
                            assembler=_get_assembler(cid),
                            contig=cid,
                            quality = QUALITY.FULL_MATCH,
                        ))
            x = _longest(_get_contigs())
            # if x.id == "CEC_1723": print(x.id, x.meta, len(x.raw_seq))
            _candidates.append(x)

        # fwd and rev hits -> possible to scaffold
        _select = lambda a, b: [(k, h) for k, hs in a.items() if k not in b for h in hs]
        ufhits, urhits = _select(fhits, rhits), _select(rhits, fhits)
        if len(ufhits)>0 and len(urhits)>0:
            _to_scaffold.append((insert_id, ufhits, urhits))

        # fwd hits only
        if len(fhits)>0:
            _candidates.append(_longest(
                Sequence(insert_id, hit.AlignSubjectToQuery(), dict(
                    assembler=_get_assembler(cid),
                    contig=cid,
                    quality = QUALITY.FORWARD_END_ONLY,
                )) for cid, hits in fhits.items() for hit in hits # this is a generator
            ))

        # rev hits only
        if len(rhits)>0:
            _candidates.append(_longest(
                Sequence(insert_id, hit.AlignSubjectToQuery(), dict(
                    assembler=_get_assembler(cid),
                    contig=cid,
                    quality = QUALITY.REVERSE_END_ONLY,
                )) for cid, hits in rhits.items() for hit in hits# this is a generator
            ))

        mapped_contigs[insert_id] = _candidates

    ######################################
    # blast contigs that mapped to fwd ends
    # against those that mapped to rev ends

    C.log.info(f"attempting to merge hits for [{len(_to_scaffold)}] inserts into scaffolds")
    MIN_SCAFFOLD_OVERLAP = 150

    @dataclass
    class ScaffoldID:
        insert_id: str
        contig_id: str

        SEP = "__"

        @classmethod
        def Parse(cls, id: str):
            toks = id.split(cls.SEP)
            if len(toks) != 2:
                C.log.error(f"failed to parse scaffolding id {id}")
                return
            return ScaffoldID(*toks)

        def Aggregate(self):
            return f"{self.insert_id}{self.SEP}{self.contig_id}"
        
        def Split(self):
            return self.insert_id, self.contig_id

    # orient contigs to the end seq and aggregate to files
    # !!! class Hits will assume this orientation
    fwd_ends_for_scaffold = open(C.out_dir.joinpath("temp.scaffolding_fwds.fa"), "w")
    rev_ends_for_scaffold = open(C.out_dir.joinpath("temp.scaffolding_revs.fa"), "w")
    fwd_scaffold_ends: dict[str, Sequence] = {}
    rev_scaffold_ends: dict[str, Sequence] = {}
    for insert_id, fhit_list, rhit_list in _to_scaffold:
        def _w(seq: RawSeq, k, file):
            file.write(f">{k}"+"\n")
            file.write(str(seq));file.write("\n")
        fwd_list = [(True, hit, fwd_ends_for_scaffold, fwd_scaffold_ends) for cid, hit in fhit_list]
        rev_list = [(False, hit, rev_ends_for_scaffold, rev_scaffold_ends) for cid, hit in rhit_list]
        for is_fwd, hit, file, d in  fwd_list+rev_list:
            cid = hit.sseqid
            con = all_contigs[cid]
            k = ScaffoldID(insert_id, cid).Aggregate()
            meta = dict(is_fwd=is_fwd)
            if con.meta is not None: meta |= con.meta
            seq_obj = Sequence(k, hit.AlignSubjectToQuery(), meta)
            _w(seq_obj.Forward(), k, file)
            d[k] = seq_obj
    fwd_ends_for_scaffold.close()
    rev_ends_for_scaffold.close()

    # make reverse end blastdb for scaffolding
    _blastdb = db_folder.joinpath("scaffolding_revs")
    _scaffold_blast_results = C.out_dir.joinpath(f"scaffolding_blast_results.tsv")
    C.shell(f"""\
    echo "# scaffold makeblastdb" >>{blast_log}
    makeblastdb \
        -dbtype nucl \
        -in {rev_ends_for_scaffold.name} \
        -out {_blastdb} >>{blast_log} 2>&1
    echo "# scaffold" >>{blast_log}
    blastn \
        -perc_identity {PIDENT_THRESH} \
        -num_threads {C.threads} \
        -query {fwd_ends_for_scaffold.name} \
        -db {_blastdb} \
        -outfmt "6 {' '.join(Hit.COLS)}" \
        -out {_scaffold_blast_results} >>{blast_log} 2>&1
    """)
    scaffolding_hits: dict[str, list[tuple[Hit, str, str]]] = {}
    try:
        _df = pd.read_csv(_scaffold_blast_results, sep="\t", header=None)
        _df.columns = Hit.COLS
        for _, row in _df.iterrows():
            h = Hit(row, fwd_scaffold_ends, rev_scaffold_ends)
            fwd_id, rev_id = [ScaffoldID.Parse(k) for k in [h.qseqid, h.sseqid]]
            if fwd_id is None or rev_id is None: continue
            fwd_insert, fwd_contig = fwd_id.Split()
            rev_insert, rev_contig = rev_id.Split()
            if fwd_insert != rev_insert: continue
            if fwd_contig == rev_contig: continue
            if h.length < MIN_SCAFFOLD_OVERLAP: continue
            scaffolding_hits[fwd_insert] = scaffolding_hits.get(fwd_insert, [])+[(h, fwd_contig, rev_contig)]
    except pd.errors.EmptyDataError:
        pass
        
    # get longest scaffolds by joining overlaps if found
    # or concatenation with "-" otherwise
    num_merged = 0
    for insert_id, fhit_list, rhit_list in _to_scaffold:
        _candidates = []
        if insert_id in scaffolding_hits:
            s = _longest(Sequence("temp", hit.Union(), dict(f=f, r=r)) for hit, f, r in scaffolding_hits[insert_id])
            meta = s.meta if s.meta is not None else {}
            _candidates.append(
                Sequence(insert_id, s.Forward(), dict(
                    assembler=f"{_get_assembler(meta.get('f'))};{_get_assembler(meta.get('r'))}",
                    contig=f"{meta.get('f')};{meta.get('r')}",
                    quality = QUALITY.FULL_SCAFFOLD,
                ))
            )
            num_merged += 1

        lfs, lrs = ([Sequence(cid, h.AlignSubjectToQuery()) for cid, h in hits] for hits in [fhit_list, rhit_list])
        for lf in lfs:
            for lr in lrs:
                _candidates.append(
                    Sequence(insert_id, lf.Forward()+GAP_CHAR+lr.ReverseCompliment(), dict(
                        assembler=f"{_get_assembler(lf.id)};{_get_assembler(lr.id)}",
                        contig=f"{lf.id};{lr.id}",
                        quality = QUALITY.GAPPED_SCAFFOLD,
                    ))
                )
        mapped_contigs[insert_id] = mapped_contigs.get(insert_id,[]) + _candidates # candidates may have been added before
    C.log.info(f"[{num_merged}] scaffolded without gaps")

    ######################################
    # resolve candidates based on quality and length
    # write resolved scaffolds and report stats

    # heurstic for choosing best-guess scaffold
    def calc_score(q: QUALITY, length: int, debug=False):
        # def _sigmoid(x): # smooth transition from favouring length to quality across the expected insert length
        #     if b == 0:
        #         return 1 if x >= expected_insert_length else 0
        #     else:
        #         return 1/(1+np.exp(-(x-expected_insert_length)/b)) # default is +- ~5kbp at 20kbp 

        def _guassian(x): # prioritize quality when in expected range
            r = expected_insert_length_range
            if r == 0: return 0
            k = 1.75 # this makes the bump drop below 0.5 at +-r
            n = (x-expected_insert_length)/(k*r) # normal distribution, 2r to cross 0.5 at +-r
            return np.exp(-(n**2))
        
        # qfrac = _sigmoid(length)
        qfrac = _guassian(length)

        qscore = q.value/QUALITY.FULL_MATCH.value
        deviation = abs(expected_insert_length-length)
        lscore = max(0, 1-(deviation/expected_insert_length))
        score = (qscore*qfrac) + (lscore*(1-qfrac))
        if debug: print(score, qfrac, qscore, lscore)
        return score
    
    num_hits = len([lst for lst in mapped_contigs.values() if len(lst)>0])
    C.log.info(f"mapped [{num_hits}] of [{len(end_seqs.insert_ids)}] paired ends")
    
    num_full_hits = 0
    mapped_file_path = C.root_workspace.joinpath("scaffolds.fna")
    report_file = C.root_workspace.joinpath("mapping_report.csv")
    full_report_file = C.out_dir.joinpath("full_mapping_report.csv")
    _rows, _full_rows = [], []
    with open(mapped_file_path, "w") as f:
        _mapped = [(k, v) for k, v in mapped_contigs.items()]
        for insert_id, _candidates in sorted(_mapped, key=lambda t: t[0]):
            if len(_candidates) == 0:
                row = (insert_id, QUALITY.NO_HITS.name.lower(), False, 0, None, None)
                _rows.append(row); _full_rows.append(row)
                continue

            def _get_meta(s: Sequence) -> tuple[QUALITY, str, str]:
                meta = s.meta if s.meta is not None else {}
                return meta["quality"], meta["assembler"], meta["contig"]

            _scored_candidates: list[tuple[float, Sequence, tuple, QUALITY, str]] = []
            for _s in _candidates:
                quality, assembler, contig_id = _get_meta(_s)
                paired = (quality == QUALITY.FULL_MATCH) or (quality == QUALITY.FULL_SCAFFOLD) or (quality == QUALITY.GAPPED_SCAFFOLD)
                seq = _s.Forward()
                orig_id = ';'.join(original_contig_id.get(_k, _k) for _k in contig_id.split(';'))
                row = (insert_id, quality.name.lower(), paired, len(seq), assembler, contig_id, orig_id)
                _full_rows.append(row)
                score = calc_score(quality, len(seq))
                # if insert_id == "case_facing_but_past":
                #     print(_s.id, score, quality, len(seq))
                #     calc_score(quality, len(seq), debug=True)
                #     print()
                _scored_candidates.append((score, _s, row, quality, assembler))
            _, best_seq, best_row, quality, assembler = sorted(_scored_candidates, key=lambda t: t[0], reverse=True)[0]
            _rows.append(best_row)
            if quality == QUALITY.FULL_MATCH: num_full_hits+=1
            f.write(f">{insert_id} {quality.name.lower()} length={len(best_seq.Forward())} {assembler.replace(';', ',')}"+"\n")
            f.write(str(best_seq.Forward()) if quality != QUALITY.REVERSE_END_ONLY else str(best_seq.ReverseCompliment())); f.write("\n")
    C.log.info(f"[{num_full_hits}] were considered a [{QUALITY.FULL_MATCH.name.lower()}]")

    _cols = "id, mapping_quality, paired, resolved_length, assemblers, contig_id, original_id".split(", ")
    report = pd.DataFrame(_rows, columns=_cols)
    report.to_csv(report_file, index=False)
    pd.DataFrame(_full_rows, columns=_cols).to_csv(full_report_file, index=False)

    Scaffolds(mapped_file_path, report_file).Save(C.expected_output)
    ClearTemp(C.out_dir)
