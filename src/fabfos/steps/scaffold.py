from __future__ import annotations
import os, sys
from dataclasses import dataclass
from typing import Any, Iterable
from enum import Enum, auto
from pathlib import Path
from typing import Callable
from Bio import SeqIO
import pandas as pd
from ..models import Scaffolds, RawContigs, EndSequences
from .common import ClearTemp, Init

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
        self.percent_query_mapped = nident/qlen*100

        # to zero index
        qlen, slen, qstart, qend, sstart, send = [x-1 for x in [qlen, slen, qstart, qend, sstart, send]]
        
        # precalculate truncation of subject by query
        if qstart > qend or sstart > send: # one flipped
            s, e = send+qstart-qlen, 0
            s = max(s, e)
        else:
            s, e = sstart-qstart+qlen, slen
            s = min(s, e)
        self._subject_trunc_loc = s, e

        # precalculate union, uses query seq for overlap
        if qstart > qend or sstart > send: # one flipped
            qa1, qb1, sa1, sb1 = 0, qend, sstart, 0
            qa2, qb2, sa2, sb2 = qlen, qstart, send, slen
        else:
            qa1, qb1, sa1, sb1 = 0, qend, send, slen
            qa2, qb2, sa2, sb2 = qlen, qstart, sstart, 0
        l1 = abs((qa1-qb1)+(sa1-sb1))
        l2 = abs((qa2-qb2)+(sa2-sb2))
        self._union_loc = (qa1, qb1, sa1, sb1) if l1>l2 else (qa2, qb2, sa2, sb2)

        self.query = query_seqs[qseqid]
        self.subject = subject_seqs[sseqid]

    def _cut(self, seq: str, start: int, end: int):
        if start > end:
            return seq[start+1:end+1:-1]
        else:
            return seq[start:end]

    def SequenceFromPairedHit(self, other: Hit):
        # orienting self as forward end hit
        fs, fe = self._subject_trunc_loc
        rs, re = other._subject_trunc_loc
        if (fs > fe) == (rs > re): # both hit are in the same direction
            return ""
        s, e = fs, rs # whoever is flipped, this should be the inner coords of the 4
        return self.query.seq+self._cut(self.subject.seq, s, e)+other.query.seq[::-1]

    def AlignSubjectToQuery(self):
        s, e = self._subject_trunc_loc
        if s == e:
            return self.query.seq
        else:
            return self.query.seq+self._cut(self.subject.seq, s, e)

    def Union(self):
        qa, qb, sa, sb = self._union_loc
        qseq = self._cut(self.query.seq, qa, qb)
        sseq = self._cut(self.subject.seq, sa, sb)
        return qseq+sseq
        
def Procedure(args):
    C = Init(args, __file__)
    raw_contigs = RawContigs.Load(C.NextArg())
    end_seqs = EndSequences.Load(C.NextArg())
    assert end_seqs.forward is not None and end_seqs.reverse is not None

    ######################################
    # load end seqs and contigs as @Sequence

    MIN_LEN = C.params.get("min_length", 1000)
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

    fwd_ends: dict[str, Sequence] = {}
    rev_ends: dict[str, Sequence] = {}
    for f, d in [
        (end_seqs.forward, fwd_ends),
        (end_seqs.reverse, rev_ends),
    ]:
        for e in SeqIO.parse(f, "fasta"):
            k = str(e.id)
            d[k] = Sequence(k, str(e.seq))

    ######################################
    # blast ends against contigs

    C.log.info(f"blasting [{len(end_seqs.insert_ids)}] ends against [{len(all_contigs)}] contigs")

    # make blastdb of contigs
    db_folder = C.out_dir.joinpath("temp.blast_dbs")
    os.makedirs(db_folder, exist_ok=True)
    _blastdb = db_folder.joinpath("all_contigs")
    blast_log = C.out_dir.joinpath("blast_log.txt")
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
            return x.seq
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
        FULL_MATCH = 1
        FULL_SCAFFOLD = 0.99
        GAPPED_SCAFFOLD = 0.7
        FORWARD_END_ONLY = 0.5
        REVERSE_END_ONLY = 0.5
        NO_HITS = 0

    C.log.info(f"mapping [{num_hits}] hits ")
    mapped_contigs: dict[str, list[Sequence]] = {}
    _to_scaffold: list[tuple[str, dict[str, Hit], dict[str, Hit]]] = []
    for i, insert_id in enumerate(end_seqs.insert_ids):
        _dictify: Callable[[dict[str, list[Hit]]], dict[str, Hit]] = lambda hit_dict: {h.sseqid:h for h in hit_dict.get(insert_id, [])}
        fhits, rhits = _dictify(fwd_hits), _dictify(rev_hits)

        # full contig match
        _candidates = []
        paired_hits = fhits.keys()&rhits.keys()
        if len(paired_hits)>0:
            def _get_contigs():
                for cid in paired_hits:
                    fwd, rev = fhits[cid], rhits[cid]
                    seq = fwd.SequenceFromPairedHit(rev)
                    yield Sequence(insert_id, seq, dict(
                        assembler=_get_assembler(cid),
                        quality = QUALITY.FULL_MATCH,
                    ))
            _candidates.append(_longest(_get_contigs()))

        # fwd and rev hits, possible to scaffold
        _select = lambda a, b: {k:h for k, h in a.items() if k not in b}
        ufhits, urhits = _select(fhits, rhits), _select(rhits, fhits)
        if len(ufhits)>0 and len(urhits)>0:
            _to_scaffold.append((insert_id, ufhits, urhits))

        # fwd hits only
        if len(fhits)>0:
            _candidates.append(_longest(
                Sequence(insert_id, hit.AlignSubjectToQuery(), dict(
                    assembler=_get_assembler(cid),
                    quality = QUALITY.FORWARD_END_ONLY,
                )) for cid, hit in fhits.items() # this is a generator
            ))

        # rev hits only
        if len(rhits)>0:
            _candidates.append(_longest(
                Sequence(insert_id, hit.AlignSubjectToQuery(), dict(
                    assembler=_get_assembler(cid),
                    quality = QUALITY.REVERSE_END_ONLY,
                )) for cid, hit in rhits.items() # this is a generator
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
    fwd_ends_for_scaffold = open(C.out_dir.joinpath("temp.scaffolding_fwds.fa"), "w")
    rev_ends_for_scaffold = open(C.out_dir.joinpath("temp.scaffolding_revs.fa"), "w")
    fwd_scaffold_ends: dict[str, Sequence] = {}
    rev_scaffold_ends: dict[str, Sequence] = {}
    for insert_id, fhits, rhits in _to_scaffold:
        def _w(seq, k, file):
            file.write(f">{k}"+"\n")
            file.write(seq);file.write("\n")
        fwd_list = [(True, hit, fwd_ends_for_scaffold, fwd_scaffold_ends) for hit in fhits.values()]
        rev_list = [(False, hit, rev_ends_for_scaffold, rev_scaffold_ends) for hit in rhits.values()]
        for is_fwd, hit, file, d in  fwd_list+rev_list:
            cid = hit.sseqid
            con = all_contigs[cid]
            k = ScaffoldID(insert_id, cid).Aggregate()
            meta = dict(is_fwd=is_fwd)
            if con.meta is not None: meta |= con.meta
            seq_obj = Sequence(k, hit.AlignSubjectToQuery(), meta)
            _w(seq_obj.seq, k, file)
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
            if h.length < MIN_SCAFFOLD_OVERLAP: continue
            scaffolding_hits[fwd_insert] = scaffolding_hits.get(fwd_insert, [])+[(h, fwd_contig, rev_contig)]
    except pd.errors.EmptyDataError:
        pass
        
    # get longest scaffolds by joining overlaps if found
    # or concatenation with "-" otherwise
    num_merged = 0
    for insert_id, fhits, rhits in _to_scaffold:
        _candidates = []
        if insert_id in scaffolding_hits:
            s = _longest(Sequence("temp", hit.Union(), dict(f=f, r=r)) for hit, f, r in scaffolding_hits[insert_id])
            meta = s.meta if s.meta is not None else {}
            _candidates.append(
                Sequence(insert_id, s.seq, dict(
                    assembler=f"{_get_assembler(meta.get('f'))};{_get_assembler(meta.get('r'))}",
                    quality = QUALITY.FULL_SCAFFOLD,
                ))
            )
            num_merged += 1

        lf, lr = (_longest(Sequence(cid, h.AlignSubjectToQuery()) for cid, h in hits.items()) for hits in [fhits, rhits])
        _candidates.append(
            Sequence(insert_id, lf.seq+"-"+lr.seq[::-1], dict(
                assembler=f"{_get_assembler(lf.id)};{_get_assembler(lr.id)}",
                quality = QUALITY.GAPPED_SCAFFOLD,
            ))
        )
        mapped_contigs[insert_id] = _candidates
    C.log.info(f"[{num_merged}] scaffolded without gaps")

    ######################################
    # write resolved contigs and report stats

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
                row = (insert_id, QUALITY.NO_HITS.name.lower(), False, 0, None)
                _rows.append(row); _full_rows.append(row)
                continue

            def _get_meta(s: Sequence) -> tuple[QUALITY, str]:
                meta = s.meta if s.meta is not None else {}
                return meta["quality"], meta["assembler"]

            _scored_candidates: list[tuple[float, Sequence, tuple, QUALITY, str]] = []
            for _s in _candidates:
                quality, assembler = _get_meta(_s)
                paired = (quality == QUALITY.FULL_MATCH) or (quality == QUALITY.FULL_SCAFFOLD) or (quality == QUALITY.GAPPED_SCAFFOLD)
                row = (insert_id, quality.name.lower(), paired, len(_s.seq), assembler)
                _full_rows.append(row)
                _scored_candidates.append((quality.value*len(_s.seq), _s, row, quality, assembler))
            _, best_seq, best_row, quality, assembler = sorted(_scored_candidates, key=lambda t: t[0], reverse=True)[0]
            _rows.append(best_row)
            if quality == QUALITY.FULL_MATCH: num_full_hits+=1
            f.write(f">{insert_id} {quality.name.lower()} length={len(best_seq.seq)} {assembler.replace(';', ',')}"+"\n")
            f.write(best_seq.seq); f.write("\n")
    C.log.info(f"[{num_full_hits}] were considered a [{QUALITY.FULL_MATCH.name.lower()}]")

    _cols = "id, mapping_quality, paired, resolved_length, assemblers".split(", ")
    report = pd.DataFrame(_rows, columns=_cols)
    report.to_csv(report_file, index=False)
    pd.DataFrame(_full_rows, columns=_cols).to_csv(full_report_file, index=False)

    Scaffolds(mapped_file_path, report_file).Save(C.expected_output)
    ClearTemp(C.out_dir)
