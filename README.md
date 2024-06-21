# FabFos: an automated pipeline for resolving inserts from pooled fosmid DNA

Tony Liu, Connor Morgan-Lang, Joe Ho, Avery Noonan, Kateryna Ievdokymenko, Zach Armstrong, Steven Hallam

## For the impatient
```bash
conda install -c txyliu \
    -c imperial-college-research-computing \
    -c bioconda -c conda-forge fabfos

# minimal
fabfos --output ./example_out \
    --interleaved --reads /.../interleaved.fastq \
    --background /.../host_background_genome.fasta \
    --pool-size 384
    
# full
fabfos --output ./example_out \
    --threads 8 \
    --verbose \
    --assembler megahit \
    --reads /.../forward_reads.fastq \
    --reverse /.../reverse_reads.fastq \
    --parity pe \
    --background /.../host_background_genome.fasta \
    --vector /.../plasmid_backbone.fasta \
    --ends /.../end_sequences.fq \
    --ends-name-regex "\w+_\d+" \
    --ends-fw-flag "FW"
```

## Overview:

Fabfos is a pipeline for resolving cloned inserts from pooled fosmid libraries using the following steps:

1. Read QC
    - Filtering of the host background with BWA [[7](#references)] and Samtools [[4](#references)]
    - Quality trimming with Trimmomatic [[3](#references)]
2. Estimate pool size
    - Reads crossing one of the two junctions between the vector backbone and insert are prepared and fed to VSEARCH [[10](#references)], which estimates the number of unique sequences as a proxy to the number of clones in the pool.
    - Alternatively, the user can provide an estimate directly.
3. Assembly
    - Choice of Megahit [[6](#references)] or Spades [[9, 8, 2](#references)]
    - _(Experimental)_ Can also use Nanopore long reads with CANU [[5](#references)]
    - Fabfos calculates assembly statistics and coverage.
4. Assessing the completeness of inserts using end sequences
    - End sequences covering both junctions between the insert and vector backbone are aligned to the assembled contigs using BLAST [[1](#references)]. Contigs with mapped end sequences for both junctions are considered complete.

## Install

Pick one of Conda, Singularity, Docker, or Manual.

### _Conda_

```bash
conda install \
    -c hallamlab \
    -c imperial-college-research-computing \
    -c bioconda \
    -c conda-forge \
    fabfos
```
_Consider using [mamba](https://mamba.readthedocs.io/en/latest/mamba-installation.html#mamba-install), a drop in multithreaded replacement for conda_

### _Singularity_

```bash
singularity pull ./fabfos.sif docker://quay.io/hallamlab/fabfos

# example run
singularity exec \
    --bind ./:/ws \
    --workdir /ws \
    ./fabfos.sif fabfos --help
```

### _Docker_

```bash
docker pull quay.io/hallamlab/fabfos

# example run
docker run -it --rm \
    -u $(id -u):$(id -g) \
    --mount type=bind,source="./",target="/ws"\
    --workdir="/ws" \
    quay.io/hallamlab/fabfos fabfos --help
```

### _Manual_

Clone this repo

Install dependencies from `.yml` file
```bash
./dev.sh --ibase
# or
conda env create --no-default-packages -n fabfos_env -f ./envs/base.yml
```

Activate the environment
```bash
conda activate fabfos_env
```
and run the source code directly
```bash
./dev.sh -r --help
# or
cd ./src
python -m fabfos --help
```

## Usage

### _Explanation of arguments_

```bash
fabfos --output ./example_out \
    --threads 8 \
    --verbose \
    --assembler megahit \
    --reads /.../forward_reads.fastq \
    --reverse /.../reverse_reads.fastq \
    --parity pe \
    --background /.../host_background_genome.fasta \
    --vector /.../plasmid_backbone.fasta \
    --ends /.../end_sequences.fq \
    --ends-name-regex "\w+_\d+" \
    --ends-fw-flag "FW"
```
Explanation of arguments:
- `--output`: The folder where Fabfos should store intermediate and output files. If it doesn't exist, Fabfos will create it.
- `--threads`: Maximum number of threads to use. Fabfos itself will only ever use 1, due to the current limitations of python
- `--verbose`: Include debug messages in printouts
- `--assembler`: The assembler and preset to use. Options are:
    - `megahit` megahit is significantly faster than spades with comparable assembly performance, however each option tends to resolve a slightly different set contigs so trying them all will produce the most complete set 
    - `spades_meta`
    - `spades_isolate`
    - `spades_sc`
- `--reads`, `--reverse`: paths to the raw reads of the pooled 
    - default is paried end fastqs
    - Fabfos should handle gzipped reads automatically
    - for single end reads `fabfos ... --reads /.../reads.fq --parity se`, leaving out `--reverse`
    - for interleaved reads `fabfos ... --reads /.../reads.fq --interleaved`, leaving out `--reverse`
- `--background`: path to genomic fasta of host background. Fabfos filters out reads that map to the host. [Example for _E. coli_ k12](https://www.ncbi.nlm.nih.gov/nuccore/U00096.2?report=fasta) 
- `--vector`: path to fasta of the vector backbone sequence. An example would be pcc1
    - used to estimate the pool size
    - can be replaced with a manual estimate, example: `--pool-size 384`, in which case Fabfos will not perform an estimate
- (optional) `--ends`: end sequences should be sequenced inward from one of the two junctions between the insert and vector backbone
    - `--ends-name-regex`: regex to pull the name of the clone from the header of the end sequence fastq. Example: `"\w+_\d+"` would get "ABC_123" from ">ABC_123_FW"
    - `--ends-fw-flag` a token that, if found within the fastq header of the end sequence, would indicate that it was sequenced from the "forward" junction. Example: `FW` would indicate that ">ABC_123_FW" is the end sequence of the "forward" juction.
    - if ommitted, assembled contigs will not be checked for completeness

### _Expected outputs_
Within the output folder specified in `--output /.../NAME` there will exist the following...
- `temp_*/` work folders for various steps and tools
- `fabfos.log` main log
- `NAME_metadata.tsv` metadata table including read stats, assembly stats, and pool size estimate
- `NAME_fosmids_*.fasta` resolved fosmids
    - `NAME_fosmids_all_contigs.fasta` all assembled contigs
    - `NAME_fosmids_both_mapped.fasta` inserts (contigs) with end sequences mapped to both ends
    `NAME_fosmids_single_mapped.fasta`inserts (contigs) with only one end mapped to an end sequence
    -  `NAME_fosmids_not_mapped` contigs with that didn't map to any end sequence
- `NAME_end_mapping.tsv` blast results from mapping end sequences to assembled contigs
- `NAME_end_mapping_failures.tsv` table of provided end sequences that didn't map to any contig

Some outputs may be ommitted if some inputs are not provided. For example, the end mapping tables will not exist if end sequences were not provided.

### _Examples_

Example 1
- Paird end reads
- `megahit` as assembler
- Fabfos estimates pool size
- End sequences
```bash
fabfos --output ./example_out \
    --assembler megahit \
    --reads /.../forward_reads.fastq \
    --reverse /.../reverse_reads.fastq \
    --background /.../host_background_genome.fasta \
    --vector /.../plasmid_backbone.fasta \
    --ends /.../end_sequences.fq \
    --ends-name-regex "\w+_\d+" \
    --ends-fw-flag "FW"
```
-----------------------------------

Example 2
- Single end reads
- `spades_isolate` as assembler
- Fabfos estimates pool size
- No end sequences
```bash
fabfos --output ./example_out \
    --assembler spades_isolate \
    --reads /.../se.fastq \
    --parity se \
    --background /.../host_background_genome.fasta \
    --vector /.../plasmid_backbone.fasta \
```
-----------------------------------

Example 3
- Interleaved reads
- `spades_meta` as assembler
- User provides pool size estimate
- End sequences
    - headers look like this `">AB001-R"` and `">AB001-L"` for left and right instead of forward and reverse
```bash
fabfos --output ./example_out \
    --assembler spades_meta \
    --interleaved \
    --reads /.../se.fastq \
    --background /.../host_background_genome.fasta \
    --pool-size 384 \
    --ends /.../end_sequences.fq \
    --ends-name-regex "\w+\d+" \
    --ends-fw-flag "L"
```
-----------------------------------


### _Command line reference_
```
usage: fabfos [-b BACKGROUND] [--vector VECTOR] [-n POOL_SIZE] [--nanopore_reads NANOPORE_READS] [--skip_correction]
              [-r READS] [-2 REVERSE] [-t {B,F}] [-i] [-p {pe,se}] [-a {spades_meta,spades_isolate,spades_sc,megahit}]
              [-o OUTPUT] [-T THREADS] [-e ENDS] [--ends-name-regex ENDS_NAME_REGEX] [--overwrite]
              [--ends-fw-flag ENDS_FW_FLAG] [-v] [--verbose] [-h]

Pipeline for filtering, assembling and organizing fosmid sequence information.

Required arguments:
  -b BACKGROUND, --background BACKGROUND
                        Path to fosmid background fasta
  --vector VECTOR       Path to vector backbone fasta, used for estimating pool size, required if --pool-size not
                        given
  -n POOL_SIZE, --pool-size POOL_SIZE
                        Estimate of number of fosmids in pool, required if --vector not given

Sequence read-specific arguments:
  -r READS, --reads READS
                        Path to the forward strand file or the interleaved paired-end file. Can be in either FastQ or
                        BAM format.
  -2 REVERSE, --reverse REVERSE
                        Path to the reverse-end read file (if applicable)
  -t {B,F}, --type {B,F}
                        Enter B if input type is BAM, F for FastQ. [ DEFAULT = 'F' ]
  -i, --interleaved     Flag indicating the reads are interleaved (i.e. forward and reverse pairs are in the same
                        file). [DEFAULT = False]
  -p {pe,se}, --parity {pe,se}
                        Specifying the sequencing chemistry used, either paired-end (pe) or single-end (se). [DEFAULT
                        = 'pe']

Nanopore-specific [development] options:
  --nanopore_reads NANOPORE_READS
                        A FASTA file containing nanopore reads to be used in assembly.
  --skip_correction     Do not perform error-correction of nanopore reads using proovread

Optional arguments:
  -a {spades_meta,spades_isolate,spades_sc,megahit}, --assembler {spades_meta,spades_isolate,spades_sc,megahit}
                        Genome assembly software to use. [DEFAULT = spades_isolate]
  -o OUTPUT, --output OUTPUT
                        path to temp. workspace [DEFAULT = /tmp]
  -T THREADS, --threads THREADS
                        The number of threads that can be used [DEFAULT = 8]
  -e ENDS, --ends ENDS  FASTA file containing fosmid ends - these will be used for alignment to fosmid contigs.
  --ends-name-regex ENDS_NAME_REGEX
                        regex for getting name of endseq., ex. "\w+_\d+" would get ABC_123 from ABC_123_FW
  --overwrite           overwrite the output directory if it already exists
  --ends-fw-flag ENDS_FW_FLAG
                        string that marks forward ends, ex. "_FW" would cause ABC_123_FW and ABC_123_RE to be assigned
                        to forward and reverse respectively.

Miscellaneous options:
  -v, --version         Print the FabFos version and exit.
  --verbose             Increase the level of verbosity in runtime log.
  -h, --help            Show this help message and exit
```

## References

Please consider citing these along with Fabfos.

1. Altschul, S.F., Gish, W., Miller, W., Myers, E.W., Lipman, D.J., 1990. **Basic local alignment search tool.** J. Mol. Biol. 215, 403–410. https://doi.org/10.1016/S0022-2836(05)80360-2

2. Bankevich, A., Nurk, S., Antipov, D., Gurevich, A.A., Dvorkin, M., Kulikov, A.S., Lesin, V.M., Nikolenko, S.I., Pham, S., Prjibelski, A.D., Pyshkin, A.V., Sirotkin, A.V., Vyahhi, N., Tesler, G., Alekseyev, M.A., Pevzner, P.A., 2012. **SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing.** J. Comput. Biol. 19, 455–477. https://doi.org/10.1089/cmb.2012.0021

3. Bolger, A.M., Lohse, M., Usadel, B., 2014. **Trimmomatic: a flexible trimmer for Illumina sequence data.** Bioinformatics 30, 2114–2120. https://doi.org/10.1093/bioinformatics/btu170

4. Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M., Li, H., 2021. **Twelve years of SAMtools and BCFtools.** GigaScience 10, giab008. https://doi.org/10.1093/gigascience/giab008

5. Koren, S., Walenz, B.P., Berlin, K., Miller, J.R., Bergman, N.H., Phillippy, A.M., 2017. **Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation.** Genome Res. 27, 722–736. https://doi.org/10.1101/gr.215087.116

6. Li, D., Liu, C.-M., Luo, R., Sadakane, K., Lam, T.-W., 2015. **MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph.** Bioinformatics 31, 1674–1676. https://doi.org/10.1093/bioinformatics/btv033

7. Li, H., Durbin, R., 2009. **Fast and accurate short read alignment with Burrows–Wheeler transform.** Bioinformatics 25, 1754–1760. https://doi.org/10.1093/bioinformatics/btp324

8. Nurk, S., Meleshko, D., Korobeynikov, A., Pevzner, P.A., 2017. **metaSPAdes: a new versatile metagenomic assembler.** Genome Res. 27, 824–834. https://doi.org/10.1101/gr.213959.116

9. Prjibelski, A., Antipov, D., Meleshko, D., Lapidus, A., Korobeynikov, A., 2020. **Using SPAdes De Novo Assembler.** Curr. Protoc. Bioinforma. 70, e102. https://doi.org/10.1002/cpbi.102

10. Rognes, T., Flouri, T., Nichols, B., Quince, C., Mahé, F., 2016. **VSEARCH: a versatile open source tool for metagenomics.** PeerJ 4, e2584. https://doi.org/10.7717/peerj.2584

