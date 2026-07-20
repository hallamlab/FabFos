from pathlib import Path
from ..models import ReadsManifest, BackgroundGenome
from .common import ClearTemp, Init, AggregateReads, Suffix

def Procedure(args):
    C = Init(args, __file__)
    reads_save, background_save = C.args
    man = ReadsManifest.Load(reads_save)
    background = BackgroundGenome.Load(background_save)
    if background.ShouldSkip():
        C.log.info(f"skipping filter")
        man.Save(C.expected_output)
        return
    C.log.info(f"filtering with {background.fasta}")
    
    count = 0
    expected = len(man.forward)+len(man.singles)
    def _filter(fwd: Path, rev: Path|None = None):
        def _is_local(p: Path):
            return p.is_relative_to(C.root_workspace) or not p.is_absolute()
        T="temp."
        if rev is None:
            sr_params = ""
            inputs = f"{fwd}"
            all_local = _is_local(fwd)
            f, r, s = None, None, Suffix(T+fwd.name, '.filtered_se')
            out_params = f">{s}"
        else:
            sr_params = "-x sr"
            inputs = f"{fwd} {rev}"
            all_local = _is_local(fwd) and _is_local(rev)
            f, r, s = Suffix(T+fwd.name, '.filtered_pe'), Suffix(T+rev.name, '.filtered_pe'), Suffix(T+fwd.name, '.filtered_se')
            out_params = f"-1 {f} -2 {r} -s {s}"
        nonlocal count; count += 1
        _log_file = C.log_file
        rm_inputs_if_local = f"&& rm {inputs}" if all_local else ""
        cmd = f"""\
            cd {C.out_dir}
            BAM=temp.bam
            minimap2 -a {sr_params} -t {C.threads} --secondary=no {background.fasta} {inputs} 2>>{_log_file} \
            | samtools sort --threads {C.threads} -o $BAM --write-index - 2>>{_log_file} \
            && samtools view -ub -f 4 -@ {C.threads} $BAM \
            | samtools fastq --verbosity 1 -N {out_params} 2>>{_log_file} \
            {rm_inputs_if_local}
        """.split("    ")
        cmd = [l for l in cmd if l != ""]
        C.log.info("\n"+f">>> batch {count} of {expected}")
        C.log.info("\n".join(l.replace("\n", "") for l in cmd))
        C.shell(" ".join(cmd))
        return [p if p is None else C.out_dir.joinpath(p) for p in [f, r, s]]

    def _count_reads(rpath: Path|None):
        if rpath is None: return 0
        #https://bioinformatics.stackexchange.com/questions/935/fast-way-to-count-number-of-reads-and-number-of-bases-in-a-fastq-file
        read_sizes = C.out_dir.joinpath("temp.readcount.txt")
        C.shell(f"""\
            cat {rpath} \
            | awk 'NR % 4 == 2' \
            | wc -cl >{read_sizes} \
        """)
        with open(read_sizes) as f:
            toks = f.readline()[:-1].strip()
            if "\t" in toks: toks = toks.split("\t")
            else: toks = [t for t in toks.split(" ") if len(t)>0]
            nr, nuc = [int(t) for t in toks]
            return nr

    fwd, rev, single = [], [], []
    starting_reads = 0
    for f, r in zip(man.forward, man.reverse):
        starting_reads += _count_reads(f)+_count_reads(r)
        ff, fr, fs = _filter(f.absolute(), r.absolute())
        fwd.append(ff)
        rev.append(fr)
        single.append(fs)

    for s in man.singles:
        starting_reads += _count_reads(s)+_count_reads(s)
        _, _, fs = _filter(s.absolute())
        single.append(fs)

    C.log.info(f"filtered {sum(len(x) for x in man.AllReads())} reads files")
    C.log.info(f"[{starting_reads}] starting reads")
    filtered = AggregateReads(fwd, rev, single, C.out_dir)
    filtered.Save(C.expected_output)
    processed_reads = sum(_count_reads(p) for g in [filtered.forward, filtered.reverse, filtered.singles] for p in g)
    C.log.info(f"{processed_reads/starting_reads*100:0.2f}% [{processed_reads}] kept")
    ClearTemp(C.out_dir)
