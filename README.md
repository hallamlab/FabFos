# FabFos: A python-driven fosmid processing pipeline

Connor Morgan-Lang, Zach Armstrong and Steven J. Hallam

## Overview:

A python pipeline for generating high-quality fosmid assemblies and interrogating fosmid pseudomolecules
 for clone identification using end-sequences.

## Download:

```
git clone git@github.com:hallamlab/FabFos.git
```

Dependencies listed below will need to be independently installed.

### Dependencies

FabFos relies on multiple other bioinformatic softwares to process the fosmids.
Please consider citing these in your published work that relies on FabFos.

For adapter and quality trimming, FabFos uses __Trimmomatic__:
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina sequence data.
Bioinformatics, 30(15), 2114–2120. http://doi.org/10.1093/bioinformatics/btu170

__MEGAHIT__ is used to assemble the QC'd reads:
Li, D., Liu, C. M., Luo, R., Sadakane, K., & Lam, T. W. (2014).
MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph.
Bioinformatics, 31(10), 1674–1676. http://doi.org/10.1093/bioinformatics/btv033

__BLAST__ is used to align the end-sequences to the fosmid contigs:
Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool.
Journal of Molecular Biology, 215(3), 403–10. http://doi.org/10.1016/S0022-2836(05)80360-2

## The Minimum information for fosmid environmental DNA (MIFFED) table

The purpose of this table if to enforce minimum information storage and enable a standard output schema for users to
easily navigate when looking for data. Since these data are not stored in a database, the user must hunt for their data
of interest. Moreover, some of the information is valuable to the processing of the fosmids and can help improve the
pipelines performance if they are provided.

The minimum information required are "Sample Name", "Project", "Human selector" (this can be synonymous with project lead),
"Vector Name", and "Screen [in silico | functional]". miffed_template.csv is provided to guide users in entering the required
and optional parameters to the pipeline.

## Outputs

The master FabFos directory needs to contain an output file called *FabFos_master_metadata.tsv*.
Additionally a FabFos_ __project__ _metadata.tsv file is created and appended to for every project within the master FabFos repository.
These files store the information provided in the MIFFED table and new information from the reads, assemblies, and fosmid-end recruitment analysis.

Final outputs include the filtered reads (**Sample** _BackboneFiltered\_R\_*.fastq),
 Trimmomatic outputs (**Sample** _paired.fa and **Sample** _unpaired.fa),
 MEGAHIT assembly (**Sample** _contigs.fasta),
 and a csv file containing Nx curve data points for the assembly (**Sample** _nx.csv).
 If a FASTA file containing fosmid ends is provided various tables are generated
 with clones mapped to contigs and a FASTA file with updated headers and sequences trimmed to suit the fosmid-end alignments.