import os, sys
from pathlib import Path
import shutil
from ..models import ReadsManifest
from .common import Init, Suffix

def Procedure(args):
    C = Init(args)
    raw_contigs_save, end_seqs_save, id_regex = C.args
    reads_save = Path(reads_save)
    reads = ReadsManifest.Load(reads_save)
    assemblers = assemblers.split(",")