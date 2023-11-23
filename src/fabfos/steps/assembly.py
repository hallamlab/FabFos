import os, sys
from pathlib import Path
import shutil
import signal
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
    
    # spades and megahit make their own logs, so console out goes to /dev/null
    C.log.info(f"performing {len(assemblers)} assemblies using [{', '.join(assemblers)}]")
    raw_contigs = {}
    _log_len = 0 
    for i, assembler_mode in enumerate(assemblers):
        _info = f"{i+1} of {len(assemblers)}: {assembler_mode}"
        _log_len = max(_log_len, len(_info))
        C.log.info(_info+" "*(_log_len-len(_info)))
        if _MOCK: continue

        mode = assembler_mode.split("_")[-1]
        out = C.out_dir.joinpath(assembler_mode)
        if out.exists(): shutil.rmtree(out)
        if "spades" in assembler_mode:
            if mode == "meta":
                klist = f"-k {' '.join(str(x) for x in range(33, 124, 20))}"
            elif mode == "isolate":
                klist = f"-k {' '.join(str(x) for x in range(67, 128, 10))}"
            else:
                klist = ""
            C.shell(f"""\
                spades.py --threads {C.threads} --{mode} \
                    {klist} \
                    -1 {fwds} -2 {revs} -s {singles} -o {out} \
                    >/dev/null 2>&1
            """)
            contigs = out.joinpath("contigs.fasta")
        else: # megahit
            if mode == "sensitive":
                preset = ""
                kmin = "--k-min 71"
                mercy = "--no-mercy"
            else:
                preset = "--presets meta-large"
                mercy = "" # yes mercy
                kmin = ""
            C.shell(f"""\
                megahit --num-cpu-threads {C.threads} \
                    {mercy} {preset} {kmin} \
                    -1 {fwds} -2 {revs} -r {singles} -o {out} \
                    >/dev/null 2>&1
            """)
            contigs = out.joinpath("final.contigs.fa")
            
        if not contigs.exists():
            C.log.warn(f"assembler [{assembler_mode}] failed")
            continue
        raw_contigs[assembler_mode] = contigs

    if _MOCK or len(raw_contigs)>0:
        RawContigs(raw_contigs).Save(C.expected_output)
    else:
        C.log.error("all assemblers failed")
