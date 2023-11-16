import os, sys
from pathlib import Path
import shutil
from ..models import ReadsManifest, RawContigs, AssemblerModes
from .common import Init, Suffix

_MOCK = False
# _MOCK = True

def Procedure(args):
    C = Init(args)
    reads_save, asm_mode_save = C.args
    reads_save = Path(reads_save)
    reads = ReadsManifest.Load(reads_save)
    assemblers = AssemblerModes.Load(asm_mode_save).modes
    if _MOCK: C.log.warn("debug mock is active, no assemblers will actually run")

    def _stringify(lst):
        return ','.join(str(p) for p in lst)
    fwds, revs, singles = [_stringify(l) for l in [reads.forward, reads.reverse, reads.single]]
    # these make their own logs, so console out goes to /dev/null
    def megahit():
        out = C.out_dir.joinpath("megahit")
        if out.exists(): shutil.rmtree(out)
        if not _MOCK:
            os.system(f"""\
                megahit --num-cpu-threads {C.threads} --no-mercy \
                    -1 {fwds} -2 {revs} -r {singles} -o {out} \
                    >/dev/null 2>&1
            """)
        return out.joinpath("final.contigs.fa")

    def spades(mode):
        name = f"spades_{mode}"
        out = C.out_dir.joinpath(name)
        if out.exists(): shutil.rmtree(out)
        if not _MOCK:
            os.system(f"""\
                spades.py --threads {C.threads} --{mode} \
                    -1 {fwds} -2 {revs} -s {singles} -o {out} \
                    >/dev/null 2>&1
            """)
        return out.joinpath("contigs.fasta")

    C.log.info(f"performing {len(assemblers)} assemblies using [{', '.join(assemblers)}]")
    raw_contigs = {}
    _log_len = 0 
    for i, assembler_mode in enumerate(assemblers):
        _info = f"{i+1} of {len(assemblers)}: {assembler_mode}"
        _log_len = max(_log_len, len(_info))
        C.log.info(_info+" "*(_log_len-len(_info)))
        if "spades" in assembler_mode:
            mode = assembler_mode.split("_")[-1]
            contigs = spades(mode)
        else:
            contigs = megahit()
        if not _MOCK and not contigs.exists():
            C.log.warn(f"assembler [{assembler_mode}] failed")
            continue
        raw_contigs[assembler_mode] = contigs

    assert len(raw_contigs)>0, "all assemblers failed"
    RawContigs(raw_contigs).Save(C.output)
