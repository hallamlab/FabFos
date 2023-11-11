import os, sys
from pathlib import Path
import shutil
from ..models import ReadsManifest
from .common import Init, Suffix

def Procedure(args):
    C = Init(args)
    assemblers, reads_save = C.args
    reads_save = Path(reads_save)
    reads = ReadsManifest.Load(reads_save)
    assemblers = assemblers.split(",")

    def _stringify(lst):
        return ','.join(str(p) for p in lst)
    fwds, revs, singles = [_stringify(l) for l in [reads.forward, reads.reverse, reads.single]]
    def megahit():
        out = C.out_dir.joinpath("megahit")
        if out.exists(): shutil.rmtree(out)
        os.system(f"""\
            megahit --num-cpu-threads {C.threads} -1 {fwds} -2 {revs} -r {singles} -o {out}
        """)

    def spades(mode):
        out = C.out_dir.joinpath(f"spades_{mode}")
        os.makedirs(out, exist_ok=True)
        os.system(f"""\
            spades.py --threads {C.threads} --{mode} -1 {fwds} -2 {revs} -s {singles} -o {out}
        """)

    C.log.info(f"performing {len(assemblers)} assemblies using [{', '.join(assemblers)}]")
    for i, assembler_mode in enumerate(assemblers):
        C.log.info(f">>>{assembler_mode} {i+1} of {len(assemblers)}")
        if "spades" in assembler_mode:
            mode = assembler_mode.split("_")[-1]
            contigs = spades(mode)
        else:
            contigs = megahit()

    C.log.info(assemblers)
