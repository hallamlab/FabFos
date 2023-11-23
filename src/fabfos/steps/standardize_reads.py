from pathlib import Path
from ..models import ReadsManifest
from ..utils import MODULE_ROOT
from .common import Init, AggregateReads, ClearTemp, Suffix

def Procedure(args):
    C = Init(args, __file__.split("/")[-1].split(".")[0])
    reads_save = Path(C.args[0])
    man = ReadsManifest.Load(reads_save)

    # unzip and make local link
    def prep(r: Path):
        unzipped_name = r.name.replace(".gz", "")
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
        fwd = [prep(p) for p in man.forward]
        rev = [prep(p) for p in man.reverse]

    if len(man.interleaved)>0:
        C.log.info(f"deinterleaving {len(man.interleaved)} interleaved reads")
        for p in man.interleaved:
            f, r = deinterleave(p) # also unzips
            fwd.append(f)
            rev.append(r)

    if len(man.single)>0:
        C.log.info(f"unzipping {len(man.single)} single end reads if zipped")
        single = [prep(p) for p in man.single]

    C.log.info(f"standardized {sum(len(x) for x in man.AllReads())} reads files")
    AggregateReads(fwd, rev, single, C.out_dir).Save(C.expected_output)
    ClearTemp(C.out_dir)
