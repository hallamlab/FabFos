# This file is part of FabFos.
# 
# FabFos is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# FabFos is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with FabFos. If not, see <https://www.gnu.org/licenses/>.

# copyright 2023 Tony Liu, Connor Morgan-Lang, Avery Noonan,
# Zach Armstrong, and Steven J. Hallam

from Bio import SeqIO
from pathlib import Path
import logging

import subprocess
import os

def TrimBackbone(ws: Path, backbone: Path, contigs: Path):
    OUT = ws.joinpath("temp_trim_vector"); os.makedirs(OUT, exist_ok=True)
    logging.info("Removing vector backbone from contigs... ")
    contigs = contigs.absolute()
    backbone = backbone.absolute()
    db = "raw_contigs"
    log = "log.txt"
    mapping = "vector_mapping.tsv"

    os.system(f"""\
        cd {OUT}
        makeblastdb -in {contigs} -dbtype nucl -out {db} 1>>{log} 2>&1
        blastn -db {db} -outfmt "6 sseqid pident length sstart send" -query {backbone} -perc_identity 85 -out ./{mapping} 1>>{log} 2>&1
    """)

    hits = {}
    with open(OUT.joinpath(mapping)) as f:
        for l in f:
            id, _, length, start, end = l.split("\t")
            hits[id] = hits.get(id, [])+[(int(start), int(end))]

    cleaned_contigs = OUT.joinpath(f"no_vector.fa")
    with open(cleaned_contigs, "w") as f:
        original = SeqIO.parse(contigs, "fasta")
        for entry in original:
            if entry.id in hits:
                # mark detected vector sequeces with x 
                _seq = list(entry.seq)
                for s, e in hits[entry.id]:
                    for i in range(s-1, e):
                        _seq[i] = "x"
                seqs = []
                _curr = []
                # remove vector and split contigs if need be
                for ch in _seq:
                    if ch == "x":
                        if len(_curr)>0:
                            seqs.append("".join(_curr))
                            _curr.clear()
                    else:
                        _curr.append(ch)
                seqs.append("".join(_curr))
                seqs = [s for s in seqs if len(s)>0]
            else:
                seqs = [str(entry.seq)]

            for i, seq in enumerate(seqs):
                f.write(f">{entry.id}_{i:02} length={len(seq)}\n")
                f.write(seq+"\n")

    logging.info("done\n")
    return str(cleaned_contigs)

def FilterMinLength(contigs: Path, min_length: int):
    logging.info(f"Filtering contigs to be at least {min_length}bp... ")
    original = SeqIO.parse(contigs, "fasta")
    filtered_contigs_path = contigs.parent.joinpath(f"contigs_{min_length}.fa")
    SeqIO.write((e for e in original if len(e.seq)>=min_length), filtered_contigs_path, "fasta")
    logging.info("done\n")
    return filtered_contigs_path

def EstimateFosmidPoolSize(
        reads: list[Path],
        backbone: Path,
        workspace: Path,
        threads: int|None=None,
    ):

    BACKBONE = "vector_backbone"
    BACKBONE_HEADER = BACKBONE
    backbone = backbone.absolute()
    reads = [r.absolute() for r in reads]
    src = Path(__file__).parent

    os.makedirs(workspace, exist_ok=True)
    original_dir = os.getcwd()
    os.chdir(workspace)
    log = ""
    def shell(cmd):
        proc = subprocess.Popen(cmd, shell=True, preexec_fn=os.setsid, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout = proc.communicate()[0].decode("utf-8")
        nonlocal log
        log += stdout+"\n"
    
    try:
        # Copies fosmid backbone senquence to new directory
        BACKBONE_SIGNATURE = None # the first 7 nucleotides of backbone
        with open(backbone) as f:
            _ = f.readline() # original header
            with open(f"{BACKBONE}.fasta", 'w') as new_f:
                new_f.write(f">{BACKBONE_HEADER}\n") # make header predictable
                l = f.readline()
                BACKBONE_SIGNATURE = l[:7]
                new_f.write(l)
                for l in f: new_f.write(l)
        shell('bwa index '+f'{BACKBONE}.fasta')

        # Aligns reads to vector backbone ***PATH TO FILES UPDATED HERE***
        read_args = " ".join(str(r) for r in reads)

        shell(f'bwa mem ./{BACKBONE}.fasta {read_args} >{BACKBONE}_aln.sam')
        shell('samtools sort '+f'{BACKBONE}_aln.sam >{BACKBONE}_aln.bam')
        shell('samtools index '+f'{BACKBONE}_aln.bam')
        shell('samtools view '+f'{BACKBONE}_aln.bam {BACKBONE_HEADER}:0-20 -o '+f'{BACKBONE}-5-020.bam')
        shell('samtools view -F 0x10 '+f'{BACKBONE}-5-020.bam -o '+'020-fwd.bam')
        shell('samtools view -f 0x10 '+f'{BACKBONE}-5-020.bam -o '+'020-rev.bam')
        shell('samtools fasta '+'020-fwd.bam > '+'020-fwd.fasta')
        shell('samtools fasta '+'020-rev.bam > '+'020-rev.fasta')

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


        uclust_cmd = [
            f"vsearch",
            "-cluster_fast", f"{BACKBONE}-5-020.fasta",
            "-uc", f"{BACKBONE}-5-020.uc",
            "-id", str(0.9),
            "-consout", f"{BACKBONE}-5-020_clusters.fasta",
            "-sizeout",
            "-centroids", f"{BACKBONE}-5-020_centroids.fasta",
        ]
        if threads is not None:
            # print(f"\nusing {threads} threads for usearch")
            uclust_cmd+= ["-threads", threads]
        shell(" ".join(uclust_cmd))

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
            # print(f"""\
            #     {total} clusters resolved and {size_ones} had only 1 member.
            #     This leaves {total-size_ones} results
            # """.replace("  ", ""))
            with open(f"results.tsv", "w") as f:
                f.write("\t".join("n_fosmids, n_fosmids_and_singletons".split(", "))+"\n")
                f.write("\t".join(str(x) for x in [estimate, total])+"\n")
            return estimate
        finally:
            final_output.close()
            hits.close()

    finally:
        with open("./log", "w") as f:
            f.write(log)
        os.chdir(original_dir)
