import os, sys
from pathlib import Path
import shutil
from ..models import ReadsManifest
from .common import ClearTemp, Init, AggregateReads, Suffix

def Procedure(args):
    C = Init(args, __file__.split("/")[-1].split(".")[0])
    reads_save = Path(C.args[0])
    man = ReadsManifest.Load(reads_save)

    adapter_folder = shutil.which("trimmomatic")
    assert adapter_folder is not None
    adapter_folder = Path(adapter_folder).parents[1].joinpath("share/trimmomatic/adapters")

    count = 0
    expected_runs = len(man.forward)+len(man.single)
    def trim(fwd: Path, rev: Path|None=None):
        T = "temp."
        TR = ".trimmed"
        inputs = fwd if rev is None else f"{fwd} {rev}"
        if rev is None:
            inputs = fwd
            outputs = C.out_dir.joinpath(Suffix(T+fwd.name, TR))
            pf, pr, singles = None, None, [outputs]
        else:
            inputs = f"{fwd} {rev}"
            outputs = [
                C.out_dir.joinpath(Suffix(T+fwd.name, f"{TR}_pe")), C.out_dir.joinpath(Suffix(T+fwd.name, f"{TR}_se")),
                C.out_dir.joinpath(Suffix(T+rev.name, f"{TR}_pe")), C.out_dir.joinpath(Suffix(T+rev.name, f"{TR}_se")),
            ]
            
            pf, pr, singles = outputs[0], outputs[2], [outputs[1], outputs[3]]
            outputs = ' '.join(str(p) for p in outputs)
        nonlocal count; count += 1
        C.log.info("\n\n"+f">>> run {count} of {expected_runs}")
        C.shell(f"""\
            trimmomatic {'SE' if rev is None else 'PE'} -threads {C.threads} \
            {inputs} {outputs} \
            ILLUMINACLIP:{adapter_folder.joinpath('TruSeq3-PE.fa')}:2:3:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
            && rm {inputs}
        """)
        return pf, pr, singles
    
    fwd, rev, single = [], [], []
    for f, r in zip(man.forward, man.reverse):
        pf, pr, s = trim(f, r)
        fwd.append(pf)
        rev.append(pr)
        single += s
    for r in man.single:
        _, _, s = trim(r)
        single += s

    C.log.info(f"trimmed {sum(len(x) for x in man.AllReads())} reads files")
    AggregateReads(fwd, rev, single, C.out_dir).Save(C.expected_output)
    ClearTemp(C.out_dir)

