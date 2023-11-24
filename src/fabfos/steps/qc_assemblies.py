import os, sys
from pathlib import Path
from Bio import SeqIO
from ..models import Assembly, ReadsManifest, RawContigs, Scaffolds, LenFilteredContigs, QCStatsForAssemblies
from .common import ClearTemp, Init

def Procedure(args):
    C = Init(args, __file__)
    reads = ReadsManifest.Load(C.NextArg())
    contigs = RawContigs.Load(C.NextArg())
    assembly_params = Assembly.Load(C.NextArg())
    _processed_contigs_save = C.NextArg()
    contig_type = C.NextArg()
    if contig_type == "scaffolds":
        processed_contigs = Scaffolds.Load(_processed_contigs_save).fasta
        pc_name = "mapped_scaffolds"
    else:
        processed_contigs = LenFilteredContigs.Load(_processed_contigs_save).kept
        pc_name = "processed_contigs"
    assemblies = [(k,contigs.contigs[k].absolute()) for k in assembly_params.modes]
    assemblies.append((pc_name, processed_contigs.absolute()))

    C.log.info(f"getting QC stats for {len(assemblies)} contig files")
    log_file = C.log_file.absolute()
    stats_manifest = {}
    for i, (asm, asm_file) in enumerate(assemblies):
        C.log.info(f"{i+1} of {len(assemblies)}: {asm}")

        mm2_inputs = []
        for f, r in zip(reads.forward, reads.reverse):
            mm2_inputs += [f, r]
        for s in reads.single:
            mm2_inputs.append(s)
        _temp_bam = f"temp.bam"
        quast_out = f"quast_{asm}"
        C.shell(f"""\
            cd {C.out_dir}
            minimap2 -a -t {C.threads} --secondary=no {asm_file} {' '.join(str(p.absolute()) for p in mm2_inputs)} 2>>{log_file} \
            | samtools sort --threads {C.threads} -o {_temp_bam} --write-index - 2>>{log_file}

            bedtools genomecov -ibam {_temp_bam} -bg >{asm}.coverage.tsv
            samtools flagstat {_temp_bam} >{asm}.stats.txt
            quast -t {C.threads} --silent -o {quast_out} {asm_file} >>{log_file} 2>&1
            cp {quast_out}/report.tsv {asm}.quast_report.tsv
            tar -cf - {quast_out} | pigz --best -p {C.threads} >{asm}.quast.tar.gz \
            && rm -r {quast_out}
        """)
        stats_manifest[asm] = asm_file

    QCStatsForAssemblies(stats_manifest).Save(C.expected_output)
    ClearTemp(C.out_dir)
