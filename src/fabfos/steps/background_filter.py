import os, sys
from pathlib import Path
import shutil
from ..models import ReadsManifest, BackgroundGenome
from .common import Init, AggregateReads, Suffix

def Procedure(args):
    C = Init(args)
    reads_save, background_save = C.args
    man = ReadsManifest.Load(reads_save)
    background = BackgroundGenome.Load(background_save)
    if background.ShouldSkip():
        C.log.info(f"skipping filter")
        man.Save(C.output)
        return
    C.log.info(background)
    
    count = 0
    expected = len(man.forward)+len(man.single)
    def _filter(fwd: Path, rev: Path|None = None):
        T="temp."
        if rev is None:
            sr_params = ""
            inputs = f"{fwd}"
            f, r, s = None, None, Suffix(T+fwd.name, '.filtered_se')
            out_params = f">{s}"
        else:
            sr_params = "-x sr"
            inputs = f"{fwd} {rev}"
            f, r, s = Suffix(T+fwd.name, '.filtered_pe'), Suffix(T+rev.name, '.filtered_pe'), Suffix(T+fwd.name, '.filtered_se')
            out_params = f"-1 {f} -2 {r} -s {s}"
        nonlocal count; count += 1
        _log_file = C.log_file.name
        cmd = f"""\
            cd {C.out_dir}
            BAM=temp.bam
            minimap2 -a {sr_params} -t {C.threads} --secondary=no {background.fasta} {inputs} 2>>{_log_file} \
            | samtools sort --threads {C.threads} -o $BAM --write-index - 2>>{_log_file} \
            && samtools view -ub -f 4 -@ {C.threads} $BAM \
            | samtools fastq --verbosity 1 -N {out_params} 2>>{_log_file} \
            && rm {inputs}
        """
        C.log.info("\n\n"+f">>> run {count} of {expected}")
        C.log.info("\n"+cmd)
        os.system(cmd)
        return [p if p is None else C.out_dir.joinpath(p) for p in [f, r, s]]

    fwd, rev, single = [], [], []
    for f, r in zip(man.forward, man.reverse):
        ff, fr, fs = _filter(f.absolute(), r.absolute())
        fwd.append(ff)
        rev.append(fr)
        single.append(fs)

    for s in man.single:
        _, _, fs = _filter(s.absolute())
        single.append(fs)

    C.log.info(f"filtered {sum(len(x) for x in man.AllReads())} reads files")
    AggregateReads(fwd, rev, single, C.out_dir).Save(C.output)
    os.system(f"rm {C.out_dir}/temp*")
