import os
from pathlib import Path
import pandas as pd
from .common import Init
from ..models import MetapathwaysArgs, Scaffolds, LenFilteredContigs
from ..models import AnnotationResults

def Procedure(args):
    C = Init(args, __file__)
    insert_type = C.NextArg()
    resolved_inserts = Scaffolds.Load(C.NextArg()).fasta if insert_type == "scaffolds" else LenFilteredContigs.Load(C.NextArg()).kept
    mpw_args = MetapathwaysArgs.Load(C.NextArg())
    assert mpw_args.reference_databases is not None
    assert mpw_args.skip_annotation == False

    mpw_out = C.root_workspace.joinpath("metapathways")
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

    _name = ".".join(resolved_inserts.name.split(".")[:-1])
    r = C.shell(f"""\
        if [ -d "{mpw_out}" ]; then rm -r {mpw_out}; fi
        metapathways run {str_args} \
        && mv {mpw_out}/{_name}/* {mpw_out} && rmdir {mpw_out}/{_name}
    """)
    if r.exit_code != 0:
        C.log.error(f"metapathways failed with exit code {r.exit_code}")
        exit(1)

    AnnotationResults(raw_results=mpw_out, contigs_used=resolved_inserts).Save(C.expected_output)
