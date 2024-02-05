from pathlib import Path
from ..models import ReadsManifest
from ..utils import MODULE_ROOT
from .common import Init, AggregateReads, ClearTemp, RemoveExt, Suffix

def Procedure(args):
    C = Init(args, __file__)
    man = ReadsManifest.Load(C.NextArg())

    # unzip and make local link
    def prep(r: Path, local_name: str):
        unzipped_name = local_name.replace(".gz", "")
        out_file = C.out_dir.joinpath("temp."+unzipped_name)
        if not r.name.endswith("gz"):
            # caution: this seems useless but 
            # filtering with trimmomatic step deletes input reads 
            # when done to reduce disk usage by intermediate files
            C.shell(f"ln -s {r} {out_file}")
        else:
            C.shell(f"""\
                pigz -p {C.threads} -dkc {r} >{out_file}
            """)
        return out_file

    def deinterleave(r: Path):
        local_file = C.out_dir.joinpath("temp."+r.name.replace(".gz", ""))
        out_f, out_r = [C.out_dir.joinpath(Suffix(local_file.name, f"_{i}")) for i in [1, 2]]
        if not r.name.endswith("gz"):
            C.shell(f"""\
                {MODULE_ROOT.joinpath("steps/deinterleave_fastq.sh")} < {r} {out_f} {out_r}
            """)
        else:
            C.shell(f"""\
                pigz -p {C.threads} -dkc {r} \
                | {MODULE_ROOT.joinpath("steps/deinterleave_fastq.sh")} {out_f} {out_r}
            """)
        return out_f, out_r

    fwd, rev, single = [], [], []
    if len(man.forward)>0:
        C.log.info(f"unzipping {len(man.forward)} paired end reads if zipped")
        for f, r in zip(man.forward, man.reverse):
            name, ext = RemoveExt(f.name)
            fwd.append(prep(f, f"{name}_1.fq"))
            rev.append(prep(r, f"{name}_2.fq"))

    if len(man.interleaved)>0:
        C.log.info(f"deinterleaving {len(man.interleaved)} interleaved reads")
        for p in man.interleaved:
            f, r = deinterleave(p) # also unzips
            fwd.append(f)
            rev.append(r)

    if len(man.singles)>0:
        C.log.info(f"unzipping {len(man.singles)} single end reads if zipped")
        single = [prep(p, p.name) for p in man.singles]

    C.log.info(f"standardized {sum(len(x) for x in man.AllReads())} read files")
    AggregateReads(fwd, rev, single, C.out_dir).Save(C.expected_output)
    ClearTemp(C.out_dir)
