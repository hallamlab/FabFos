import os, sys
from pathlib import Path
import shutil
from ..models import ReadsManifest, RawContigs, AssemblerModes
from .common import Init, Suffix

def Procedure(args):
    C = Init(args)
    reads_save, asm_mode_save = C.args
    reads_save = Path(reads_save)
    reads = ReadsManifest.Load(reads_save)
    assemblers = AssemblerModes.Load(asm_mode_save)

    def _stringify(lst):
        return ','.join(str(p) for p in lst)
    fwds, revs, singles = [_stringify(l) for l in [reads.forward, reads.reverse, reads.single]]
    def megahit():
        out = C.out_dir.joinpath("megahit")
        if out.exists(): shutil.rmtree(out)
        os.system(f"""\
            megahit --num-cpu-threads {C.threads} -1 {fwds} -2 {revs} -r {singles} -o {out}
        """)
        return out.joinpath("final.contigs.fa")

    def spades(mode):
        out = C.out_dir.joinpath(f"spades_{mode}")
        if out.exists(): shutil.rmtree(out)
        os.makedirs(out, exist_ok=True)
        os.system(f"""\
            spades.py --threads {C.threads} --{mode} -1 {fwds} -2 {revs} -s {singles} -o {out}
        """)
        return out.joinpath("contigs.fasta")

    C.log.info(f"performing {len(assemblers)} assemblies using [{', '.join(assemblers)}]")
    raw_contigs = {}
    for i, assembler_mode in enumerate(assemblers):
        C.log.info(f">>>{assembler_mode} {i+1} of {len(assemblers)}")
        if "spades" in assembler_mode:
            mode = assembler_mode.split("_")[-1]
            contigs = spades(mode)
        else:
            contigs = megahit()
        if not contigs.exists():
            C.log.warn(f"assembler [{assembler_mode}] failed")
            continue
        raw_contigs[assembler_mode] = contigs

    assert len(raw_contigs)>0, "all assemblers failed"
    RawContigs(raw_contigs).Save(C.output)
