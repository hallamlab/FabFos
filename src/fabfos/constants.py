from pathlib import Path

READS_FOLDER = Path("temp_reads")
ORIGINAL_READS = READS_FOLDER.joinpath("original_reads.json")

ASM_FOLDER = Path("temp_assembly")
ASM_MODES = ASM_FOLDER.joinpath("assemblers.json")

ENDSEQ_FOLDER = Path("temp_endmap")
ENDSEQS = ENDSEQ_FOLDER.joinpath("endseqs.json")

SKIP_FILTER = "skip_filter"
