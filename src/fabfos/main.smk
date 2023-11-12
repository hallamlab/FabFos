import os
from pathlib import Path
from fabfos.models import ReadsManifest, BackgroundGenome, AssemblerModes, EndSequences
from fabfos.models import RawContigs

# the actual commands to call each tool are found in src/fabfos/steps/*.py

SRC = Path(config["src"])
LOGS = Path(config["log"])
THREADS = int(config["threads"])

rule all:
    input: f"{ASSEMBLY}"

rule standardize_reads:
    input: f"{ORIGINAL_READS}"
    output: f"{STD_READS}"
    params:
        src=SRC
    threads: THREADS
    shell:
        """\
        PYTHONPATH={params.src}:$PYTHONPATH
        python -m fabfos api --step standardize_reads --args {threads} {output} {input}
        """

rule quality_trim:
    input: f"{STD_READS}"
    output: f"{TRIMMED_READS}"
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
        f"{TRIMMED_READS}",
        f"{config['background']}"
    output: f"{FILTERED_READS}"
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
        f"{FILTERED_READS}",
        f"{ASM_MODES}"
    output: f"{ASSEMBLY}"
    params:
        src=SRC,
    threads: THREADS
    shell:
        """\
        PYTHONPATH={params.src}:$PYTHONPATH
        python -m fabfos api --step assembly --args {threads} {output} {input}
        """

# rule end_map:
#     input: f"{ASSEMBLY}"
#     params:
#         src=SRC,
#         assemblers=config['assemblers']
#     threads: THREADS