import os
from pathlib import Path
from fabfos.models import ReadsManifest, BackgroundGenome, Assembly, EndSequences
from fabfos.models import RawContigs, Scaffolds, LenFilteredContigs

# the actual commands to call each tool are found in src/fabfos/steps/*.py

SRC = Path(config["src"])
LOGS = Path(config["log"])
THREADS = int(config["threads"])

GIVEN_ENDSEQS = config["endseqs"]

    # input: f"{RawContigs.MANIFEST}"
rule all:
    input: f"{Scaffolds.MANIFEST if GIVEN_ENDSEQS else LenFilteredContigs.MANIFEST}"

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
        f"{Assembly.ARG_FILE}"
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
# Scaffolds or LenFilteredContigs

# if given end seqs
rule scaffold:
    input:
        f"{RawContigs.MANIFEST}",
        f"{EndSequences.ARG_FILE}"
    output: f"{Scaffolds.MANIFEST}"
    params:
        src=SRC,
    threads: THREADS
    shell:
        """\
        PYTHONPATH={params.src}:$PYTHONPATH
        python -m fabfos api --step scaffold --args {threads} {output} {input}
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