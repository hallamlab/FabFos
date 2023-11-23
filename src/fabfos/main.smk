import os
from pathlib import Path
from fabfos.models import ReadsManifest, BackgroundGenome, AssemblerModes, EndSequences
from fabfos.models import RawContigs, EndMappedContigs, LenFilteredContigs

# the actual commands to call each tool are found in src/fabfos/steps/*.py

SRC = Path(config["src"])
LOGS = Path(config["log"])
THREADS = int(config["threads"])

GIVEN_ENDSEQS = config["endseqs"]

    # input: f"{RawContigs.MANIFEST}"
rule all:
    input: f"{EndMappedContigs.MANIFEST if GIVEN_ENDSEQS else LenFilteredContigs.MANIFEST}"

rule standardize_reads:
    input: f"{ReadsManifest.ARG_FILE}"
    output: f"{ReadsManifest.STD}"
    params:
        src=SRC
    threads: THREADS
    shell:
        """\
        PYTHONPATH={params.src}:$PYTHONPATH
        python -m fabfos api --step standardize_reads --args {threads} {output} {input}
        """

rule quality_trim:
    input: f"{ReadsManifest.STD}"
    output: f"{ReadsManifest.TRIM}"
    params:
        src=SRC,
    threads: THREADS
    shell:
        """\
        PYTHONPATH={params.src}:$PYTHONPATH
        python -m fabfos api --step quality_trim --args {threads} {output} {input}
        """

rule filter_background:
    input:
        f"{ReadsManifest.TRIM}",
        f"{BackgroundGenome.ARG_FILE}"
    output: f"{ReadsManifest.FILTER}"
    params:
        src=SRC,
    threads: THREADS
    shell:
        """\
        PYTHONPATH={params.src}:$PYTHONPATH
        python -m fabfos api --step background_filter --args {threads} {output} {input}
        """

rule assembly:
    input:
        f"{ReadsManifest.FILTER}",
        f"{AssemblerModes.ARG_FILE}"
    output: f"{RawContigs.MANIFEST}"
    params:
        src=SRC,
    threads: THREADS
    shell:
        """\
        PYTHONPATH={params.src}:$PYTHONPATH
        python -m fabfos api --step assembly --args {threads} {output} {input}
        """

# ---------------------------------
# conditional branching by requesting either
# EndMappedContigs or LenFilteredContigs

# if given end seqs
rule select_contigs:
    input:
        f"{RawContigs.MANIFEST}",
        f"{EndSequences.ARG_FILE}"
    output: f"{EndMappedContigs.MANIFEST}"
    params:
        src=SRC,
    threads: THREADS
    shell:
        """\
        PYTHONPATH={params.src}:$PYTHONPATH
        python -m fabfos api --step select_contigs --args {threads} {output} {input}
        """

# else
rule filter_contigs_by_length:
    input: f"{RawContigs.MANIFEST}"
    output: f"{LenFilteredContigs.MANIFEST}"
    params:
        src=SRC,
    threads: THREADS
    shell:
        """\
        PYTHONPATH={params.src}:$PYTHONPATH
        python -m fabfos api --step filter_contigs_by_length --args {threads} {output} {input}
        """
# end if
# ---------------------------------