import os, sys
from pathlib import Path
import shutil
from ..models import ReadsManifest, RawContigs, Assembly
from .common import ClearTemp, Init, Suffix

_MOCK = False
# _MOCK = True

def Procedure(args):
    C = Init(args, __file__)
    reads_save, asm_mode_save = C.args
    reads_save = Path(reads_save)
    reads = ReadsManifest.Load(reads_save)
    asm_meta = Assembly.Load(asm_mode_save)
    assemblers = asm_meta.modes
    given_contigs = asm_meta.given
    if _MOCK: C.log.warn("debug mock is active, no assemblers will actually run")

    def _stringify(lst):
        return ','.join(str(p) for p in lst)
    fwds, revs, singles = [_stringify(l) for l in [reads.forward, reads.reverse, reads.single]]
    
    # spades and megahit make their own logs, so console out goes to /dev/null
    if len(assemblers)>0: C.log.info(f"performing {len(assemblers)} assemblies using [{', '.join(assemblers)}]")
    if len(given_contigs)>0: C.log.info(f"registering {len(given_contigs)} given assemblies")
    raw_contigs = {}
    contigs_dir = asm_meta.CONTIG_DIR
    os.makedirs(contigs_dir, exist_ok=True)

    _expected_len = len(assemblers) + len(given_contigs)
    # assembled_contigs = []
    for i, assembler_mode in enumerate(list(assemblers)+list(given_contigs.keys())):
        if _MOCK: continue

        mode = assembler_mode.split("_")[-1]
        asm_out = C.out_dir.joinpath("temp."+assembler_mode)
        expected_out = contigs_dir.joinpath(f"{assembler_mode}.fna")
        if expected_out.exists():
            C.log.info(f"{i+1} of {_expected_len}: existing file [{assembler_mode}] registered")
            raw_contigs[assembler_mode] = expected_out
            continue
        else:
            C.log.info(f"{i+1} of {_expected_len}: [{assembler_mode}]")
            if asm_out.exists(): shutil.rmtree(asm_out)

        if assembler_mode not in Assembly.CHOICES:
            C.log.error(f"[{assembler_mode}] contigs not present and not a known assembler mode, skipping")
            continue

        if "spades" in assembler_mode:
            if mode == "meta":
                klist = f"-k {' '.join(str(x) for x in range(33, 124, 20))}"
            elif mode == "isolate":
                klist = f"-k {' '.join(str(x) for x in range(67, 128, 10))}"
            else:
                klist = ""
            r = C.shell(f"""\
                spades.py --threads {C.threads} --{mode} \
                    {klist} \
                    -1 {fwds} -2 {revs} -s {singles} -o {asm_out} \
                    >/dev/null 2>&1 \
                && cp {asm_out}/contigs.fasta {expected_out} \
                && cp {asm_out}/spades.log {C.root_workspace}/logs/assembly.{assembler_mode}.log
            """)
        else: # megahit
            if mode == "sensitive":
                preset = ""
                kmin = "--k-min 71"
                mercy = "--no-mercy"
            else:
                preset = "--presets meta-large"
                mercy = "" # yes mercy
                kmin = ""
            r = C.shell(f"""\
                megahit --num-cpu-threads {C.threads} \
                    {mercy} {preset} {kmin} \
                    -1 {fwds} -2 {revs} -r {singles} -o {asm_out} \
                    >/dev/null 2>&1 \
                && cp {asm_out}/final.contigs.fa {expected_out} \
                && cp {asm_out}/log {C.root_workspace}/logs/assembly.{assembler_mode}.log
            """)
            if r.killed: return
            
        if not expected_out.exists():
            C.log.warn(f"assembler [{assembler_mode}] failed")
            continue
        raw_contigs[assembler_mode] = expected_out
        # assembled_contigs.append(assembler_mode)

    if _MOCK or len(raw_contigs)>0:
        RawContigs(raw_contigs).Save(C.expected_output)
        ClearTemp(C.out_dir)
    else:
        C.log.error("all assemblers failed")
