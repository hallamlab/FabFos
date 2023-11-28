import os
from pathlib import Path
from Bio import SeqIO
import pandas as pd
from .common import Init
from ..models import ReadsManifest, VectorBackbone
from ..models import PoolSizeEstimate
from ..process_management import Shell

def Procedure(args):
    C = Init(args, __file__)
    reads = ReadsManifest.Load(C.NextArg())
    vector = VectorBackbone.Load(C.NextArg())
    assert vector.fasta is not None, "vector backbone fasta not given"
    assert len(reads.forward) == 1, "reads were not aggregated"

    fwd, rev = reads.forward[0], reads.reverse[0]

    BACKBONE_SIGNATURE_SIZE = 7
    BACKBONE = "vector_backbone"
    BACKBONE_HEADER = BACKBONE
    backbone = vector.fasta
    workspace = C.out_dir

    os.makedirs(workspace, exist_ok=True)
    original_dir = os.getcwd()
    os.chdir(workspace)

    def _shell(cmd: str):
        def _log(x: str):
            with open(C.log_file, "a") as f:
                f.write(x);
        # r = Shell(cmd, lambda x: C.log.info(x), lambda x: C.log.info(f"ERR: {x}"))
        r = Shell(cmd, _log, lambda x: _log(f"STD_ERR: {x}"))
        if r.killed:
            C.log.error("killed")
            exit(1)
    try:
        # Copies fosmid backbone senquence to new directory
        BACKBONE_SIGNATURE = None # the first 7 nucleotides of backbone
        with open(backbone) as f:
            _ = f.readline() # original header
            with open(f"{BACKBONE}.fasta", 'w') as new_f:
                new_f.write(f">{BACKBONE_HEADER}\n") # make header predictable
                l = f.readline()
                BACKBONE_SIGNATURE = l[:BACKBONE_SIGNATURE_SIZE]
                new_f.write(l)
                for l in f: new_f.write(l)

        # Aligns reads to vector backbone
        # assumes fwd ends in _1 and rev ends in _2
        _shell(f"""\
            BAM=temp.bam
            minimap2 -a -x sr -t {C.threads} --secondary=no ./{BACKBONE}.fasta {fwd} {rev} \
            && samtools view -ub -F 4 -@ {C.threads} $BAM \
            | samtools fastq --verbosity 1 -N -1 020-fwd.fasta -2 020-rev.fasta
        """)

        fwd = SeqIO.parse(str('020-fwd.fasta'), "fasta")
        recs = []
        for f in fwd:
            recs.append((True, f))

        rev = SeqIO.parse(str('020-rev.fasta'), "fasta")
        for r in rev:
            rev_r = r.seq.reverse_complement()
            r.seq = rev_r
            recs.append((False, r))

        kept_indicies = []
        final_recs = []
        hit_lens = []
        CUT = 100
        for i, (is_fwd, rec) in enumerate(recs):
            if BACKBONE_SIGNATURE not in rec.seq: continue
            new_seq = rec.seq[:rec.seq.index(BACKBONE_SIGNATURE)]
            hit_lens.append(len(new_seq))
            if len(new_seq) >= CUT:
                trim_seq = new_seq[len(new_seq)-CUT:len(new_seq)]
                final_recs.append((is_fwd, trim_seq))
                kept_indicies.append(i)
        assert len(final_recs)>0, "no hits to vector backbone"
        with open(f'{BACKBONE}-5-020.fasta', 'w') as f:
            for i, (is_fwd, seq) in enumerate(final_recs):
                f.write(f">{i}\n{seq}\n")

        with open(f"hits_info.tsv", "w") as f:
            f.write("\t".join("hit_count, hits_kept, hit_lengths".split(", "))+"\n")
            hl_str = ";".join(str(l) for l in hit_lens)
            f.write("\t".join(str(x) for x in [len(hit_lens), len(final_recs), hl_str])+"\n")


        _shell(f"""
            vsearch -sizeout -id 0.9 -threads {C.threads} \
            -cluster_fast {BACKBONE}-5-020.fasta \
            -uc {BACKBONE}-5-020.uc \
            -consout {BACKBONE}-5-020_clusters.fasta \
            -centroids {BACKBONE}-5-020_centroids.fasta
        """)

        final_output =  open(f'original.fasta', 'w')
        hits =          open(f'hits.fasta', 'w')
        try:
            centroids = SeqIO.parse(str(f'{BACKBONE}-5-020_centroids.fasta'), "fasta")
            size_ones, estimate = 0, 0
            for c in centroids:
                id, size = c.id.split(";")
                id = int(id)
                rev, original = recs[kept_indicies[id]]
                header = f">{original.id};{'forward' if not rev else 'reverse_compliment'};{size}"
                final_output.write(f"{header}\n{original.seq}\n")
                hits.write(f"{header}\n{c.seq}\n")

                estimate += 1
                if size == "size=1":
                    size_ones += 1

            total = estimate + size_ones

            os.chdir(original_dir)
            C.log.info(f"estimated pool size: {estimate}")
            pd.DataFrame([(estimate, total)], columns=["estimated_pool_size", "estimated_pool_size_with_singletons"])\
                .to_csv(C.root_workspace.joinpath("pool_size_estimate.csv"), index=False)
            C.log.info(C.expected_output)
            PoolSizeEstimate(size=estimate, size_with_singletons=total).Save(C.expected_output)

        finally:
            final_output.close()
            hits.close()

    finally:
        os.chdir(original_dir)
