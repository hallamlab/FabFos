from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable
import json

@dataclass
class ReadsManifest:
    forward: list[Path]
    reverse: list[Path]
    interleaved: list[Path]
    single: list[Path]

    def AllReads(self):
        return [self.forward, self.reverse, self.interleaved, self.single]

    def Verify(self, on_error: Callable):
        all_reads = self.AllReads()
        if sum(len(x) for x in all_reads) == 0:
            on_error("no reads given")
        for p in [p for g in all_reads for p in g]:
            if p.exists(): continue
            on_error(f"[{p}] does not exist")
        if len(self.forward) != len(self.reverse):
            on_error("number of forward and reverse reads don't match")

    def Save(self, path: Path):
        with open(path, "w") as j:
            json.dump({k:[str(v) for v in l] for k, l in self.__dict__.items()}, j, indent=4)

    @classmethod
    def Load(cls, path: Path):
        with open(path) as j:
            raw = {k: [Path(v) for v in l] for k, l in json.load(j).items()}
            return ReadsManifest(**raw)
