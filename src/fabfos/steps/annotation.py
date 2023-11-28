import os
from pathlib import Path
import pandas as pd
from .common import Init
from ..models import MetapathwaysArgs, Scaffolds, LenFilteredContigs
from ..models import MetapathwaysResults
from ..process_management import Shell

def Procedure(args):
    C = Init(args, __file__)
    insert_type = C.NextArg()
    resolved_inserts = Scaffolds.Load(C.NextArg()).fasta if insert_type == "scaffolds" else LenFilteredContigs.Load(C.NextArg()).kept
    mpw_args = MetapathwaysArgs.Load(C.NextArg())
    assert mpw_args.reference_databases is not None
    assert mpw_args.skip_annotation == False

    def _shell(cmd: str):
        def _log(x: str):
            with open(C.log_file, "a") as f:
                f.write(x);
        r = Shell(cmd, _log, lambda x: _log(f"STD_ERR: {x}"))
        if r.killed:
            C.log.error("killed")
            exit(1)
        return r

    mpw_out = C.root_workspace.joinpath("annotation")
    args = dict(
        input_file = resolved_inserts,
        refdb_dir = mpw_args.reference_databases,
        COMPUTE_TPM = "skip",
        output_dir = mpw_out,
    )
    functional_dbs_path = mpw_args.reference_databases.joinpath("functional/formatted")
    if functional_dbs_path.exists():
        dbs = [f.replace("-names.txt", "") for f in os.listdir(functional_dbs_path) if "names.txt" in f]
        args |= dict(annotation_dbs = ' '.join(dbs))

    manual_overrides = mpw_args.additional_args
    BL = ["output_dir", "o"]
    for k in BL:
        if k in manual_overrides:
            C.log.error(f"arg [{k}] is not allowed to be set manually and wil be ignored")
            del manual_overrides[k]
    args |= manual_overrides

    str_args = " ".join([f"--{k} {v}" for k,v in args.items()])

    r = _shell(f"""\
        metapathways run {str_args}
    """)
    if r.exit_code == 0:
        MetapathwaysResults(mpw_out).Save(C.expected_output)
