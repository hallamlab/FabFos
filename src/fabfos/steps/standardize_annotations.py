import os
from pathlib import Path
import json
from dataclasses import dataclass, field
from typing import Iterable
import pandas as pd
from pytest import skip
from .common import Init
from ..models import AnnotationResults
from ..models import StandardizedAnnotations
from ..utils import Batchify

def Procedure(args):
    C = Init(args, __file__)
    mpw_results = AnnotationResults.Load(C.NextArg())
    
    mpw_out = mpw_results.raw_results
    ann_out = C.root_workspace.joinpath("annotations")
    C.shell(f"""\
        if [ -d "{ann_out}" ]; then rm -r {ann_out}; fi
        mkdir {ann_out}
    """)

    _name = mpw_results.contigs_used.name
    NAME = ".".join(_name.split(".")[:-1]) # name of input fasta file ("contigs", or "scaffolds" depending on previous steps)

    contig_mapping_path = mpw_out.joinpath(f"preprocessed/{NAME}.mapping.txt")
    mpw2original = {}
    for _, row in pd.read_csv(contig_mapping_path, sep="\t", header=None, names=["contig", "original_name", "length"]).iterrows():
        mpw2original[row.contig[len(f"{NAME}-"):]] = row.original_name

    # -------------------------------------------------------------------
    # parse open reading frame fastas 
    # and filter based on metapathways qc

    @dataclass
    class Orf:
        contig: str
        orf: str
        aa: str
        nt: str
    
    def _parse_fasta(path: Path):
        entries = {}
        with open(mpw_out.joinpath(path)) as f:
            _contig, _orf, _seq = "", "", []
            def _submit():
                nonlocal _contig, _orf, _seq
                if len(_seq)==0: return
                entries[(_contig, _orf)] = "".join(_seq)
                _contig, _orf = "", ""
                _seq.clear()

            for l in f:
                l = l[:-1]
                if l.startswith(">"):
                    _submit()
                    _contig, _orf = l[1:].split("-")[1:]
                    _orf = _orf[len("G"):]
                    _seq = []
                else:
                    _seq.append(l)
            _submit()
        return entries
    aa_orfs = _parse_fasta(mpw_out.joinpath(f"orf_prediction/{NAME}.qced.faa"))
    nt_orfs = _parse_fasta(mpw_out.joinpath(f"orf_prediction/{NAME}.fna"))
    all_orfs, mapping = {}, {}
    max_orf_count = 0
    last_contig, i = None, 1
    for k, seq in aa_orfs.items():
        c, o = k
        orig_contig = mpw2original[c]
        if last_contig == orig_contig:
            i += 1
        else:
            last_contig, i = orig_contig, 1
        max_orf_count = max(max_orf_count, i)
        new_orf = str(i) # leading zeros added shortly below
        mapping[(c, o)] = orig_contig, new_orf
        all_orfs[(orig_contig, new_orf)] = Orf(orig_contig, new_orf, seq, nt_orfs[k])
    orf_places = len(str(max_orf_count))
    def _prefix(o: str):
        return f"{int(o):0{orf_places}}"
    def _prefix_orf(orf: Orf):
        orf.orf = _prefix(orf.orf)
        return orf
    all_orfs = {(c, _prefix(o)):_prefix_orf(orf) for (c, o), orf in all_orfs.items()}
    mapping = {(c, o):(oc, _prefix(no)) for (c, o), (oc, no) in mapping.items()}
    pd.DataFrame([(oc, no, c, o) for (c, o), (oc, no) in mapping.items()], columns=["contig", "orf", "metapathways_contig", "metapathways_orf"])\
        .to_csv(ann_out.joinpath("metapathways_mapping.csv"), index=False)
    
    # -------------------------------------------------------------------
    # get CDS as template for gff3s and sort
    # write ordered orf fastas

    @dataclass
    class GffCdsTemplate:
        contig: str
        start: int
        end: int
        is_plus: bool
        phase: int

        def MakeLine(self, id: str, name: str|None=None, score: float|None=None, attributes: dict=dict()):
            # https://useast.ensembl.org/info/website/upload/gff3.html
            # http://gmod.org/wiki/GFF3
            # http://www.sequenceontology.org/browser/obob.cgi

            final_attributes = dict(
                ID=id,
            )|attributes
            if name is not None:
                final_attributes|=dict(
                    Name=name,
                )
            clean = lambda x: x.replace("=", ":").replace(";", "|")
            return "\t".join([
                self.contig, ".", "CDS",                                # seqid, source, type
                str(self.start), str(self.end),                         # start, end
                str(score) if score is not None else ".",               # score
                "+" if self.is_plus else "-",                           # strand
                str(self.phase),                                        # phase
                ";".join([f"{clean(k)}={clean(v)}" for k,v in final_attributes.items()])    # attributes
            ])
    
    templates: dict[tuple[str, str], GffCdsTemplate] = {}
    with open(mpw_out.joinpath(f"orf_prediction/{NAME}.cds.gff")) as f:
        for l in f:
            if l.startswith("#"): continue
            _contig, _, _, _start, _end, _, strand, _phase, _attributes = l.split("\t")
            start, end, phase = [int(x) for x in [_start, _end, _phase]]
            is_plus = strand == "+"
            attributes = dict([x.split("=") for x in _attributes.split(";")[:-1]])
            contig = _contig[len(f"{NAME}-"):]
            orf = attributes["ID"].split("_")[1]
            templates[(contig, orf)] = GffCdsTemplate(mpw2original[contig], start, end, is_plus, phase)

    gffs = {}
    def write_gff(name: str, entries: Iterable[str]):
        path = ann_out.joinpath(f"{name}.gff")
        with open(path, "w") as f:
            f.write(f"# GFF3 - {name}\n")
            for line in entries:
                f.write(line); f.write("\n")
        gffs[name] = path

    ordered_templates = []
    for k, _template in templates.items():
        if k not in mapping: continue
        contig, orf = k
        orig_contig, new_orf = mapping[(contig, orf)]
        ordered_templates.append((f"{orig_contig} {min(_template.start, _template.end):032}", orig_contig, new_orf, _template))
    ordered_templates.sort(key=lambda x: x[0])

    templates.clear()
    _entries = []
    _fwds, _revs = [], []
    for _, orig_contig, new_orf, _template in ordered_templates:
        templates[(orig_contig, new_orf)] = _template
        line = _template.MakeLine(id=f"{orig_contig}_{new_orf}")
        if _template.is_plus:
            _fwds.append(line)
        else:
            _revs.append(line)
        _entries.append(line)
    write_gff("open_reading_frames", _entries)
    write_gff("open_reading_frames_fwd", _fwds)
    write_gff("open_reading_frames_rev", _revs)

    # write orf sequences with correct order and headers
    faa_path = ann_out.joinpath("open_reading_frames.faa")
    faa = open(faa_path, "w")
    fna_path = ann_out.joinpath("open_reading_frames.fna")
    fna = open(fna_path, "w")
    for _, orig_contig, new_orf, _ in ordered_templates:
        _orf = all_orfs[(orig_contig, new_orf)]
        faa.write(f">{orig_contig}_{new_orf}\n")
        fna.write(f">{orig_contig}_{new_orf}\n")
        faa.write(_orf.aa+"\n")
        fna.write(_orf.nt+"\n")
        # for seq in Batchify(_orf.aa, 80): faa.write(str(seq)+"\n")
        # for seq in Batchify(_orf.nt, 80): fna.write(str(seq)+"\n")
    faa.close()
    fna.close()

    # -------------------------------------------------------------------
    # parse functionl annotation hits as gff3

    @dataclass
    class Hit:
        ref_id: str
        ref_desc: str
        bsr: float

    blast_results_path = mpw_out.joinpath("blast_results")
    for file_name in [f for f in os.listdir(blast_results_path) if f.endswith(".parsed.txt")]:
        ref_name = file_name.split(".")[1]
        hits_dict: dict[tuple[str, str], list[Hit]] = {}
        for _, row in pd.read_csv(blast_results_path.joinpath(file_name), sep="\t").iterrows():
            q, ref_id, ref_desc, bsr = [row[k] for k in "#query target product bsr".split(" ")]
            contig, orf = q.split("-")
            contig, orf = mapping[(contig, orf[len("G"):])]
            hits_dict[(contig, orf)] = hits_dict.get(contig, []) + [Hit(ref_id, ref_desc, bsr)]
        
        _entries = []
        for (contig, orf), hits in hits_dict.items():
            if (contig, orf) not in templates:
                C.log.error(f"{ref_name} had hit that was not in predicted orfs") # should not be possible
                continue
            
            best_hit = hits[0]
            for candidate in hits:
                if candidate.bsr > best_hit.bsr:
                    best_hit = candidate
            
            _template = templates[(contig, orf)]
            attributes = dict(
                ref_id=best_hit.ref_id,
            )
            if len(hits)>1:
                attributes["other_hits"]=f"[{','.join([h.ref_id for h in hits])}]"
            line = _template.MakeLine(
                id=f"{contig}_{orf}",
                name=best_hit.ref_desc,
                score=min(float(best_hit.bsr), 1),
                attributes=attributes
            )
            _entries.append((f"{contig} {min(_template.start, _template.end):032}", line))
        _entries.sort(key=lambda x: x[0])
        write_gff(ref_name, [l for _, l in _entries])

    # -------------------------------------------------------------------
    # rRNA gff3

    rrna_gff_path = mpw_out.joinpath(f"orf_prediction/{NAME}.rRNA.gff")
    rrna_fna_path = mpw_out.joinpath(f"orf_prediction/{NAME}.rRNA.fna")
    rrna_results = os.listdir(mpw_out.joinpath("results/rRNA"))
    _filter = lambda lst, term: [f for f in lst if term in f.lower()]
    ssu_results, lsu_results = [_filter(rrna_results, term)[0] for term in ["ssu", "lsu"]]

    # -------------------------------------------------------------------
    # tRNA gff3

    trna_gff_path = mpw_out.joinpath(f"orf_prediction/{NAME}.tRNA.gff")
    trna_fna_path = mpw_out.joinpath(f"orf_prediction/{NAME}.tRNA.fasta")
    trna_results_path = mpw_out.joinpath(f"results/tRNA/{NAME}.tRNA.results.txt")
    trna_results = pd.read_csv(
        trna_results_path, sep="\t", skiprows=3, header=None,
        names="contig trna_number start end aa codon intron_start intron_end score note".split(" ")
    )

    StandardizedAnnotations(
        gffs = {k:p for k,p in gffs.items()},
        fna = fna_path,
        faa = faa_path,
    ).Save(C.expected_output)
