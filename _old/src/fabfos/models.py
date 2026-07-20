import os
from dataclasses import dataclass, field
from pathlib import Path
from re import L
from typing import Callable, Any
from Bio import SeqIO
import json
import pandas as pd

from .utils import regex

# CREATED_FOLDERS:set[Path] = set()
# # this automatically registers the first folder of static paths in saveable models
# class AutoRegisterPaths(type):
#     def __new__(cls, name, bases, dct):
#         x = super().__new__(cls, name, bases, dct)
#         for v in x.__dict__.values():
#             if isinstance(v, Path) and not v.is_absolute():
#                 d = v.parents[-2] if len(v.parents) > 1 else v
#                 CREATED_FOLDERS.add(d)
#         return x

# class Saveable(metaclass=AutoRegisterPaths):
class Saveable:
    def Save(self, path: str|Path):
        path = Path(path)
        if not path.parent.exists(): os.makedirs(path.parent)

        def _can_save(k, v):
            if k.upper() == k: return False
            if callable(v): return False
            if isinstance(k, str) and k[0] == "_": return False
            return True

        def _stringyfy(v):
            if isinstance(v, list):
                return [str(x) for x in v]
            elif isinstance(v, dict):
                return {k:_stringyfy(x) for k, x in v.items()}
            else:
                return str(v)
        with open(path, "w") as j:
            json.dump(_stringyfy({k:v for k, v in self.__dict__.items() if _can_save(k, v)}), j, indent=4)

@dataclass
class ReadsManifest(Saveable):
    forward: list[Path]
    reverse: list[Path]
    interleaved: list[Path]
    singles: list[Path]

    ARG_FILE =  Path("internals/temp_reads/original_reads.json")
    STD =       Path("internals/temp_reads/std_reads.json")
    TRIM =      Path("internals/temp_trim/trimmed.json")
    FILTER =    Path("internals/temp_filter/filtered.json")

    def AllReads(self):
        return [self.forward, self.reverse, self.interleaved, self.singles]

    @classmethod
    def Parse(cls, args, out_dir: Path, on_error: Callable):
        _listify = lambda arg: [Path(p).absolute() for p in arg]
        model = cls(
            forward=_listify(args.forward),
            reverse=_listify(args.reverse),
            interleaved=_listify(args.interleaved),
            singles=_listify(args.single),
        )

        # verify
        all_reads = model.AllReads()
        if sum(len(x) for x in all_reads) > 0:
            for p in [p for g in all_reads for p in g]:
                if p.exists(): continue
                on_error(f"[{p}] does not exist")
            if len(model.forward) != len(model.reverse):
                on_error("number of forward and reverse reads don't match")
        
        model.Save(out_dir.joinpath(cls.ARG_FILE))
        return model

    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = {k: [Path(v) for v in l] for k, l in json.load(j).items()}
            return cls(**raw)

@dataclass
class BackgroundGenome(Saveable):
    fasta: Path|str

    ARG_FILE = ReadsManifest.FILTER.parent.joinpath("background.json")
    SKIP = "SKIP"

    def ShouldSkip(self):
        return self.fasta == self.SKIP

    @classmethod
    def Parse(cls, args, out_dir: Path, on_error: Callable):
        if args.background is not None:
            bg = Path(args.background).absolute()
            if not bg.exists(): on_error(f"[{bg}] does not exist")
        else:
            bg = cls.SKIP
        model = cls(bg)
        model.Save(out_dir.joinpath(cls.ARG_FILE))
        return model

    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = {k:Path(v) if v != cls.SKIP else str(v) for k, v in json.load(j).items()}
            return cls(**raw)
        
@dataclass
class Assembly(Saveable):
    modes: list[str]
    given: dict[str, Path]
    CHOICES = "megahit_sensitive, megahit_default, megahit_meta, spades_isolate, spades_meta, spades_sc".lower().split(", ")

    ARG_FILE = Path("internals/temp_assembly/assemblies.json")
    CONTIG_DIR = Path("intermediate_contigs")

    @classmethod
    def Parse(cls, args, out_dir: Path, on_error: Callable):
        picked_assemblers = []
        _seen = set()
        _given = {}
        for a in args.assemblies:
            if a in _seen: continue
            if a not in cls.CHOICES:
                asm_path = Path(a)
                if not asm_path.exists():
                    on_error(f"[{a}] doesn't exist and is not an assembler mode")
                else:
                    k = asm_path.name
                    _given[k] = _given.get(k, [])+[asm_path]
            else:
                a = str(a).lower()
                if a not in _seen: picked_assemblers.append(a)
            _seen.add(a)

        contig_folder = out_dir.joinpath(cls.CONTIG_DIR)
        os.makedirs(contig_folder, exist_ok=True)
        given_contigs = {}
        def _remove_extension(fname: str):
            name_toks = fname.split(".")
            return ".".join(name_toks[:-1]) if len(name_toks)>1 else name_toks[0]
        def _make_local(apath: Path, name: str):
            local_path = contig_folder.joinpath(f"{name}.fna")
            if local_path.exists(): os.unlink(local_path)
            os.link(apath, local_path)
            given_contigs[name] = apath
        for asm_lst in _given.values():
            if len(asm_lst) == 1:
                a = asm_lst[0]
                _make_local(a, f"given_{_remove_extension(a.name)}")
            else:
                for i, apath in enumerate(asm_lst):
                    _make_local(apath, f"{_remove_extension(apath.name)}_{i+1:04}")
        
        model = cls(picked_assemblers, given_contigs)
        model.Save(out_dir.joinpath(cls.ARG_FILE))
        return model

    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = json.load(j)
            return cls(**raw)

@dataclass
class EndSequences(Saveable):
    given: bool
    forward: Path|None
    reverse: Path|None
    insert_ids: list[str] # as in fosmid inserts

    MANIFEST_COLUMNS = ["file", "sequence_id", "name", "is_forward"]
    ARG_FILE = Path("internals/temp_scaffold/endseqs.json")
    SKIP = "SKIP"

    @classmethod
    def Parse(cls, args, out_dir: Path, on_error: Callable):
        NULL_MODEL = cls(False, None, None, [])
        save_path = out_dir.joinpath(cls.ARG_FILE)
        NULL_MODEL.Save(save_path)
        if args.ends_manifest is None:
            return NULL_MODEL

        try:
            ends_man = Path(args.ends_manifest)
        except ValueError:
            on_error("given ends manifest not a valid file path")
            return NULL_MODEL
        if not ends_man.exists():
            on_error(f"given path to ends manifest doesn't point to a file")
            return NULL_MODEL
        
        file_type = ends_man.suffix
        if file_type == ".csv":
            dfman = pd.read_csv(ends_man)
        elif file_type == ".tsv":
            dfman = pd.read_csv(ends_man, sep="\t")
        elif file_type in {".xlsx", ".xls"}:
            dfman = pd.read_excel(ends_man)
        else:
            on_error(f"ends manifest file type [{file_type}] not recognized")
            return NULL_MODEL

        columns = dfman.columns
        for c in cls.MANIFEST_COLUMNS:
            if c not in columns:
                on_error(f"ends manifest missing column [{c}]")
                return NULL_MODEL
            
        positive = ["yes", "y", "true", 1]
        negative = ["no", "n", "false", 0]
        pos_neg = positive+negative
        for _, r in dfman.iterrows():
            is_forward = r["is_forward"]
            if isinstance(is_forward, str): is_forward = is_forward.lower()
            if is_forward not in pos_neg:
                on_error(f"ends manifest is_forward column must be in the form: yes/no, y/n, true/false, 1/0, got [{is_forward}]")
                return NULL_MODEL
            
        WS = cls.ARG_FILE.parent
        if not WS.exists(): os.makedirs(WS)

        seqs = {}
        _err = False        
        for fpath in dfman["file"].unique():
            try:
                file_path = Path(fpath)
                if not file_path.exists():
                    on_error(f"ends file [{fpath}] doesn't exist")
                    _err = True; continue
            except ValueError:
                on_error(f"ends file [{fpath}] is not valid")
                _err = True; continue
            ext = file_path.suffix
            if ext in {".fasta", ".fa", ".fna"}:
                ext = "fasta"
            elif ext in {".fastq", ".fq"}:
                ext = "fastq"
            elif ext in {".ab1"}:
                ext = "ab1"
            else:
                on_error(f"ends file [{file_path}] has unrecognized extension [{ext}]")
                _err = True; continue
            for e in SeqIO.parse(file_path, ext):
                seqs[fpath, str(e.id)] = e
        if _err: return NULL_MODEL

        allf_p, allr_p = [out_dir.joinpath(cls.ARG_FILE.parent).joinpath(f).absolute() for f in ["endf.fa", "endr.fa"]]
        if not allf_p.parent.exists(): os.makedirs(allf_p.parent)
        allf, allr = [open(p, "w") for p in [allf_p, allr_p]]

        fids = []
        rids = []
        for _, r in dfman.iterrows():
            is_forward = r["is_forward"]
            if isinstance(is_forward, str): is_forward = is_forward.lower()
            out = fids if is_forward in positive else rids
            file_path = r["file"]
            header = str(r["sequence_id"])
            out.append((str(r["name"]), (Path(file_path), seqs[(file_path, header)])))

        # _no_pair = lambda p, x: f"{p} [{x}] has no matching pair"
        _dup = lambda p, x: f"{p} [{x}] is duplicate"
        try:
            unmatched = 0
            all_ids = set()
            for this, other, out in [
                (fids, rids, allf),
                (rids, fids, allr),
            ]:  
                _seen = set()
                other = set(id for id, _ in other)
                for id, (p, e) in this:
                    all_ids.add(id)
                    if id in _seen: on_error(_dup(p, id))
                    if id not in other: unmatched += 1
                    _seen.add(id)
                    out.write(f">{id}"+"\n")
                    out.write(str(e.seq))
                    out.write("\n")
            if unmatched > 0: print(f"WARNING: [{unmatched}] ends had no matching pair")
        finally:
            allf.close()
            allr.close()
        
        model = cls(True, allf_p, allr_p, list(all_ids))
        model.Save(save_path)
        return model
    
    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = json.load(j)
            _pathify = lambda x: Path(x) if x != "None" else None
            return cls(
                given = raw.get("given", "False").title() == "True",
                forward = _pathify(raw.get("forward", "None")),
                reverse = _pathify(raw.get("reverse", "None")),
                insert_ids = raw.get("insert_ids", []),
            )

@dataclass
class VectorBackbone(Saveable):
    vector_backbone_given: bool
    fasta: Path|None

    ARG_FILE = Path("internals/temp_estimate_pool_size/manifest.json")
    SKIP = "SKIP"

    @classmethod
    def Parse(cls, args, out_dir: Path, on_error: Callable):
        raw_vec_path = args.vector
        if raw_vec_path is not None:
            raw_vec_path = Path(raw_vec_path).absolute()
            if not raw_vec_path.exists():
                on_error(f"vector backbone file [{raw_vec_path}] doesn't exist")
            
        model = cls(
            vector_backbone_given=raw_vec_path is not None,
            fasta=raw_vec_path,
        )
        model.Save(out_dir.joinpath(cls.ARG_FILE))
        return model

    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = json.load(j)
            _pathify = lambda x: Path(x) if x != "None" else None
            return cls(
                vector_backbone_given = raw.get("vector_backbone_given", "False").title() == "True",
                fasta = _pathify(raw.get("fasta", "None")),
            )
        
@dataclass
class MetapathwaysArgs(Saveable):
    skip_annotation: bool
    reference_databases: Path|None
    additional_args: dict[str, str]

    ARG_FILE = Path("internals/temp_annotation_raw/metapathways_args.json")

    @classmethod
    def Parse(cls, args, out_dir: Path, on_error: Callable):
        if args.reference_databases is None:
            model = cls(
                skip_annotation=True,
                reference_databases=None,
                additional_args={},
            )
        else:
            ref_path = Path(args.reference_databases).absolute()
            if not ref_path.exists():
                on_error(f"metapathways reference databases at [{ref_path}] don't exist")
            model = cls(
                skip_annotation=False,
                reference_databases=ref_path,
                additional_args={k:v for k, v in [a.split("=") if "=" in a else (a, "") for a in args.metapathways]},
            )
        model.Save(out_dir.joinpath(cls.ARG_FILE))
        return model
    
    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = json.load(j)
            ref_path = raw.get("reference_databases", "None")
            if ref_path == "None":
                ref_path = None
            else:
                ref_path = Path(ref_path)
            return cls(**dict(
                skip_annotation = raw.get("skip_annotation", "False").title() == "True",
                reference_databases = ref_path,
                additional_args = raw.get("additional_args", {}),
            ))

#################################
# internal (not parsed from args)
#################################

@dataclass
class RawContigs(Saveable):
    contigs: dict[str, Path]

    MANIFEST = Path("internals/temp_assembly/contigs.json")

    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = json.load(j)
            _contigs = raw.get("contigs")
            if not isinstance(_contigs, dict): _contigs = {}
            return cls({k: Path(v) for k, v in _contigs.items()})

@dataclass
class Scaffolds(Saveable):
    fasta: Path
    report: Path

    MANIFEST = EndSequences.ARG_FILE.parent.joinpath("end_mapped_contigs.json")

    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = {k: Path(v) for k, v in json.load(j).items()}
            return cls(**raw)
        
@dataclass
class LenFilteredContigs(Saveable):
    kept: Path
    discarded: Path

    MANIFEST = EndSequences.ARG_FILE.parent.joinpath("len_filtered_contigs.json")

    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = {k: Path(v) for k, v in json.load(j).items()}
            return cls(**raw)

@dataclass
class QCStatsForAssemblies(Saveable):
    qcd_contigs: dict[str, Path]

    MANIFEST = Path("qc_assemblies/manifest.json")

    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = {str(k):Path(v) for k, v in json.load(j).items()}
            return cls(raw)

@dataclass
class QCStatsForReads(Saveable):
    qcd_reads: list[Path]

    MANIFEST = Path("qc_reads/manifest.json")

    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = [Path(v) for v in json.load(j)]
            return cls(raw)

@dataclass
class PoolSizeEstimate(Saveable):
    size: int
    size_with_singletons: int

    MANIFEST = Path("internals/temp_estimate_pool_size/size.json")

    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = json.load(j)
            return cls(**raw)

@dataclass
class AnnotationResults(Saveable):
    raw_results: Path
    contigs_used: Path

    MANIFEST = Path("internals/temp_annotation_raw/manifest.json")

    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = json.load(j)
            _pathify = lambda k:  Path(raw[k])
            return cls(**dict(
                raw_results = _pathify("raw_results"),
                contigs_used = _pathify("contigs_used"),
            ))
        
@dataclass
class StandardizedAnnotations(Saveable):
    gffs: dict[str, Path]
    fna: Path
    faa: Path

    MANIFEST = Path("internals/temp_annotation_std/manifest.json")

    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = json.load(j)
            return cls(
                gffs = {k:Path(p) for k, p in raw["gffs"].items()},
                fna = Path(raw["fna"]),
                faa = Path(raw["faa"]),
            )
