from ..models import ReadsManifest, QCStatsForReads
from .common import ClearTemp, Init

def Procedure(args):
    C = Init(args, __file__)
    reads = ReadsManifest.Load(C.NextArg())

    all_reads = [r for g in reads.AllReads() for r in g]
    C.shell(f"""\
        fastqc --noextract -o {C.out_dir} {" ".join([str(r) for r in all_reads])} 2>/dev/null
    """)

    QCStatsForReads(all_reads).Save(C.expected_output)
