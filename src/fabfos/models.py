import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Any
from Bio import SeqIO
import json

from .utils import regex

class Saveable:
    def Save(self, path: Path):
        if not path.parent.exists(): os.makedirs(path.parent)

        def _can_save(k, v):
            if k.upper() != k: return False
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

    ARG_FILE =  Path("temp_reads/original_reads.json")
    STD =       Path("temp_reads/std_reads.json")
    TRIM =      Path("temp_trim/trimmed.json")
    FILTER =    Path("temp_filter/filtered.json")

    def AllReads(self):
        return [self.forward, self.reverse, self.interleaved, self.single]

    @classmethod
    def Parse(cls, args, on_error: Callable):
        _listify = lambda arg: [Path(p).absolute() for p in arg]
        model = cls(
            forward=_listify(args.forward),
            reverse=_listify(args.reverse),
            interleaved=_listify(args.interleaved),
            single=_listify(args.single),
        )

        # verify
        all_reads = model.AllReads()
        if sum(len(x) for x in all_reads) == 0:
            on_error("no reads given")
        for p in [p for g in all_reads for p in g]:
            if p.exists(): continue
            on_error(f"[{p}] does not exist")
        if len(model.forward) != len(model.reverse):
            on_error("number of forward and reverse reads don't match")
        return model

    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = {k: [Path(v) for v in l] for k, l in json.load(j).items()}
            return cls(**raw)

@dataclass
class BackgroundGenome(Saveable):
    background: Path|str

    ARG_FILE = ReadsManifest.FILTER.parent.joinpath("background.json")
    SKIP = "SKIP"

    @classmethod
    def Parse(cls, args, on_error: Callable):
        if args.background is not None:
            bg = Path(args.background).absolute()
            if not bg.exists(): on_error(f"[{bg}] does not exist")
        else:
            bg = cls.SKIP
        return cls(bg)

    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = {k:Path(v) if v != cls.SKIP else str(v) for k, v in json.load(j).items()}
            return cls(**raw)
        
@dataclass
class AssemblerModes(Saveable):
    modes: list[str]
    CHOICES = "megahit, spades_meta, spades_isolate, spades_sc,".lower().split(", ")

    ARG_FILE = Path("temp_assembly/assemblers.json")

    def __len__(self):
        return len(self.modes)
    
    def __iter__(self):
        return self.modes.__iter__()

    @classmethod
    def Parse(cls, args, on_error: Callable):
        picked_assemblers = []
        _seen = set()
        for a in args.assemblers:
            a = str(a).lower()
            if a in _seen: continue
            if a not in cls.CHOICES:
                on_error(f"[{a}] is not one of {cls.CHOICES}")
            picked_assemblers.append(a)
            _seen.add(a)
        return cls(picked_assemblers)

    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = json.load(j)
            return cls(raw)

@dataclass
class EndSequences(Saveable):
    given: bool
    forward: Path|None
    reverse: Path|None
    id_regex: str|None

    ARG_FILE = Path("temp_endmap/endseqs.json")
    SKIP = "SKIP"
    DEFAULT_REGEX = r'\w+'

    def _get_id(self, s):
        return next(regex(self.id_regex, s))

    @classmethod
    def Parse(cls, args, on_error: Callable):
        _pathify = lambda p: Path(p) if p is not None else p
        model = cls(
            given = args.endf is not None,
            forward = _pathify(args.endf),
            reverse = _pathify(args.endr),
            id_regex = args.id_regex if args.id_regex is not None else cls.DEFAULT_REGEX,
        )
        paths = [model.given, model.reverse]
        if len([v for v in paths if v is None])==1:
            on_error(f"both --endf and --endr must be given or omitted together")

        if not model.given: return model # skip id match check

        assert model.forward is not None and model.reverse is not None
        for p in paths:
            if p.exists(): continue
            on_error(f"[{p}] doesn't exist")

        def _get(p):
            for e in SeqIO.parse(p, "fasta"):
                yield model._get_id(e.id)
        
        _no_pair = lambda dir, x: f"{dir} [{x}] has no matching pair"
        _dup = lambda dir, x: f"{dir} [{x}] is duplicate"
        
        _seen = set()
        fids = list(_get(model.forward))
        rids = list(_get(model.reverse))
        sids = set(rids)
        for x in fids:
            if x in _seen: on_error(_dup("endf", x))
            if x not in sids: on_error(_no_pair("endf", x))
            _seen.add(x)

        _seen.clear()
        sids = set(fids)
        for x in rids:
            if x in _seen: on_error(_dup("endr", x))
            if x not in sids: on_error(_no_pair("endr", x))
            _seen.add(x)
        return model
    
    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw: dict = {k:Path(v) if "id" not in k else v for k, v in json.load(j).items()}
            return cls(**raw)

#################################
# internal (not parsed from args)
#################################

@dataclass
class RawContigs(Saveable):
    contigs: dict[str, Path]

    SAVE = Path("temp_assembly/contigs.json")

    @classmethod
    def Load(cls, path):
        with open(path) as j:
            raw = {k: Path(v) for k, v in json.load(j).items()}
            return cls(raw)
