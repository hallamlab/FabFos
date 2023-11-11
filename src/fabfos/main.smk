import os
from pathlib import Path
from fabfos.constants import READS_FOLDER, ORIGINAL_READS

# the actual commands to call each tool are found in src/fabfos/steps/*.py

SRC = Path(config["src"])
LOGS = Path(config["log"])
THREADS = int(config["threads"])

STD_READS = READS_FOLDER.joinpath("reads.json")
TRIMMED_READS = Path("temp_trim/trimmed.json")
FILTERED_READS = Path("temp_filter/filtered.json")
ASSEMBLY = Path("temp_assembly/assembly.json")

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
    input: f"{FILTERED_READS}",
    output: f"{ASSEMBLY}"
    params:
        src=SRC,
        assemblers=config['assemblers']
    threads: THREADS
    shell:
        """\
        PYTHONPATH={params.src}:$PYTHONPATH
        python -m fabfos api --step assembly --args {threads} {output} {params.assemblers} {input}
        """
