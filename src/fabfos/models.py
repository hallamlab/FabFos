import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Any
from Bio import SeqIO
import json

from .utils import regex

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
    single: list[Path]

    ARG_FILE =  Path("internals/temp_reads/original_reads.json")
    STD =       Path("internals/temp_reads/std_reads.json")
    TRIM =      Path("internals/temp_trim/trimmed.json")
    FILTER =    Path("internals/temp_filter/filtered.json")

    def AllReads(self):
        return [self.forward, self.reverse, self.interleaved, self.single]

    @classmethod
    def Parse(cls, args, out_dir: Path, on_error: Callable):
        _listify = lambda arg: [Path(p).absolute() for p in arg]
        model = cls(
            forward=_listify(args.forward),
            reverse=_listify(args.reverse),
            interleaved=_listify(args.interleaved),
            single=_listify(args.single),
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
    CHOICES = "megahit_sensitive, megahit_meta, spades_meta, spades_isolate, spades_sc,".lower().split(", ")

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
                _make_local(a, _remove_extension(a.name))
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

    ARG_FILE = Path("internals/temp_contigs/endseqs.json")
    SKIP = "SKIP"
    DEFAULT_REGEX = r'\w+'

    @classmethod
    def Parse(cls, args, out_dir: Path, on_error: Callable):
        _pathify = lambda l: [Path(p).absolute() for p in l] if l is not None else None
        forwards = _pathify(args.endf)
        reverses = _pathify(args.endr)
        id_regex = args.end_regex if args.end_regex is not None else cls.DEFAULT_REGEX
        save_path = out_dir.joinpath(cls.ARG_FILE)

        if len([v for v in [forwards, reverses] if v is None])==1:
            on_error(f"both --endf and --endr must be given or omitted together")

        if forwards is None or reverses is None: # skip id match check
            model = cls(False, None, None, [])
            model.Save(save_path)
            return model

        # ensure ids are in pairs and unique + aggregate to 2 files
        for e, p in [(e, p) for e, l in [("endf", forwards), ("endr", reverses)] for p in l]:
            if p.exists(): continue
            on_error(f"--{e} file doesn't exist [{p}]")
        
        _no_pair = lambda p, x: f"{p} [{x}] has no matching pair"
        _dup = lambda p, x: f"{p} [{x}] is duplicate"
        
        allf_p, allr_p = [out_dir.joinpath(cls.ARG_FILE.parent).joinpath(f).absolute() for f in ["endf.fa", "endr.fa"]]
        if not allf_p.parent.exists(): os.makedirs(allf_p.parent)
        allf, allr = [open(p, "w") for p in [allf_p, allr_p]]

        def _get(lst):
            for p in lst:
                if not p.exists(): continue
                for e in SeqIO.parse(p, "fasta"):
                    yield next(regex(id_regex, e.id)), (p, e)
        fids = list(_get(forwards))
        rids = list(_get(reverses))
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
                if id not in other: on_error(_no_pair(p, id))
                _seen.add(id)
                out.write(f">{id}"+"\n")
                out.write(str(e.seq))
                out.write("\n")
            out.close()
        
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
