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

try:
    import argparse
    import sys
    import glob
    import os
    import re
    import subprocess
    import shutil
    import logging
    # import time
    from Bio import SeqIO
    from pathlib import Path
    from packaging import version
    from .addons import EstimateFosmidPoolSize, TrimBackbone, FilterMinLength
    from .external_qc import Fastqc, AssemblyStats
    
except (ImportWarning, ModuleNotFoundError):
    import traceback
    sys.stderr.write("Could not load some user defined module functions")
    sys.stderr.write(str(traceback.print_exc(10)))
    sys.exit(3)

HERE = Path(os.path.realpath(__file__)).parent
with open(HERE.joinpath("version.txt")) as f:
    VERSION = f.read()

class FabFos:
    def __init__(self, path):
        self.workspace = path

        # Add new information:
        self.os = os_type()
        if sys.version_info > (2, 9):
            self.py_version = 3
        else:
            self.py_version = 2

        # this assumes a conda environment
        self.adaptor_trim = f"{Path(sys.executable).parents[1]}/share/trimmomatic/adapters/"
        self.metadata_header = "#Sample Name (LLLLL-PP-WWW)\tProject\tHuman selector\tVector Name\t" \
                               "Screen [in silico | functional]\tSelection criteria\tNumber of fosmids\t" \
                               "Sequencing submission date (YYYY-MM-DD)\tGlycerol plate name\tSequencing center\t" \
                               "Sequencing type\tRead length\tInstrument\tDate of FabFos analysis\tNumber of reads\t" \
                               "% off-target reads\tNumber of trimmed reads\tAssembler version\tNumber of Contigs\t" \
                               "N50\tN90\tContigs > 27kbp\tContigs > 50kbp\t" \
                               "# single-pair\t# multi-pair\t# single-orphan\t# multi-orphan\t# Unidentifiable\t" \
                               "# Fosmid ends(X2 = #reads)\t# aligned fosmid ends\t# unaligned fosmid ends\t# Failed ends\n"
        return

    def furnish(self) -> None:
        os.makedirs(self.workspace, exist_ok=True)
        return

class FosmidEnds:
    """
    Class for holding count data for fosmid ends
    """
    def __init__(self):
        self.fasta_path = ""
        self.ends = None
        self.failed = set()
        self.all_clones = set()
        self.aligned = set()
        self.unaligned = set()

        self.single_pair = 0
        self.multi_pair = 0
        self.single_orphan = 0
        self.multi_orphan = 0
        self.unassigned = 0
        self.total_ends = 0
        self.num_missing = 0
        self.num_stunted = 0
        self.num_total = 0
        self.num_failed = 0
        self.num_aligned = 0
        self.num_unaligned = 0

    def load_ends(self):
        self.ends = read_fasta(self.fasta_path)  # type: dict
        self.get_fosmid_ends_stats()
        return

    def get_fosmid_ends_stats(self, min_length=100) -> None:
        if not self.ends:
            logging.error("Fosmid ends FASTA file has not been loaded - FosmidsEnds.ends is empty.\n")
            sys.exit(13)

        for name, seq in self.ends.items():
            if name[0] == '>':
                clone_name = get_fosmid_end_name(name[1:])
            else:
                clone_name = get_fosmid_end_name(name)
            self.all_clones.add(clone_name)
            if len(seq) < min_length:
                self.failed.add(clone_name)

        self.num_total = len(self.all_clones)
        self.num_failed = len(self.failed)
        return


class FosmidClone:
    def __init__(self, contig, sequence):
        self.contig = contig
        self.clone = ""
        self.evidence = ""
        self.ends = ""
        self.strand = ""
        self.sequence = sequence
        self.seq_length = 0

    def extract_sequence(self, start, end):
        if end < start:
            self.sequence = reverse_complement(self.sequence)
            self.sequence = self.sequence[end-1:start-1]
        else:
            self.sequence = self.sequence[start-1:end-1]
        self.seq_length = len(self.sequence)

    def get_info(self):
        info_string = "contig = " + str(self.contig) + \
                      "\tclone = " + str(self.clone) + \
                      "\tevidence = " + str(self.evidence) + \
                      "\tends = " + str(self.ends) + \
                      "\tstrand = " + str(self.strand) + \
                      "\tlength = " + str(self.seq_length)
        return info_string


class ReadStats:
    def __init__(self):
        self.num_raw_reads = 0  # Number of reads read from FASTQ file(s)
        self.num_on_target = 0  # Number of reads removed during filtering
        self.num_filtered_reads = 0  # Number of reads remaining after filtering
        self.num_reads_trimmed = 0  # Number of reads remaining after quality filtering
        self.num_reads_assembled = 0  # Number of reads provided to the assembler

    def calc_on_target(self) -> None:
        self.num_on_target = self.num_raw_reads - self.num_filtered_reads
        return

    def percent_on_target(self) -> float:
        return round(float(self.num_on_target*100) / self.num_raw_reads, ndigits=4)


class Sample:
    def __init__(self, sample_id):
        # General sample information
        self.id = sample_id
        self.read_stats = ReadStats()
        self.output_dir = ""
        self.assembled_fosmids = ""
        self.assembler = ""
        self.assembly_mode = ""
        self.num_fosmids_estimate = 0
        self.parity = ""
        self.interleaved = False  # Boolean indicating whether the reads in FASTQ are interleaved (True) or not (False)
        self.forward_reads = None
        self.reverse_reads = None
        self.pe_trimmed = []
        self.se_trimmed = []

        # Control flow
        self.overwrite = False
        self.assemble = True
        self.map_ends = True
        self.exclude = False

        # Nanopore-specific
        self.nanopore: str|bool = ""
        self.error_correction = True

    def retrieve_separate_pe_fq(self):
        forward = None
        reverse = None
        for fastq_file in self.pe_trimmed:
            if re.search(r'pe.1.fq$', fastq_file):
                forward = fastq_file
            elif re.search(r'pe.2.fq$', fastq_file):
                reverse = fastq_file
            else:
                logging.error("Unrecognized FASTQ file name: {}!\n".format(fastq_file))
                sys.exit(3)
        if not forward or not reverse:
            logging.error("Unable to find the paired-end FASTQ files for assembling" + self.id + ".\n")
            sys.exit(3)
        return forward, reverse

    def gather_reads(self, reads_file: str, reverse_file: str, parity: str,
                     executables: dict, output_dir: str, file_type: str) -> None:
        """
        Finds all the FASTQ files for a sample, filling the Sample.forward_reads and Sample.reverse_reads attributes
        If the file_type == 'B' indicating BAM format, these files are also converted to FASTQs

        :param reads_file: Path to the forward-oriented or interleaved FASTQ files or BAMs
        :param reverse: Path the the reverse-oriented FASTQ files
        :param parity: String indicating the sequencing chemistry library for the library [pe (default)|se]
        :param executables: Dictionary containing paths to executables that is indexed by the executable name
        :param output_dir: Used only if the input files are BAMs - path to write the new FASTQ files
        :param file_type: String indicating whether the input reads are provided in FASTQ (F) or BAM (B) format
        :return: None
        """

        folder = str(Path(self.output_dir).joinpath("temp_reads"))
        os.makedirs(folder, exist_ok=True)
        if file_type == "B":  # If the input are BAM files, convert them to FASTQ and split them here
            fastq_list = bam2fastq(reads_file, self.id, folder, executables["samtools"], parity)
            raw_reads = {}
            for fastq in fastq_list:
                if re.search("1.fastq", fastq):
                    raw_reads["forward"] = fastq
                else:
                    raw_reads["reverse"] = fastq
        else:
            raw_reads = {"forward": reads_file}
            if reverse_file is not None: raw_reads["reverse"] = reverse_file
        if self.interleaved:
            raw_reads["forward"], raw_reads["reverse"] = deinterleave_fastq(raw_reads["forward"], folder)
            self.interleaved = False

        self.forward_reads = raw_reads["forward"]
        self.reverse_reads = raw_reads.get("reverse")
        self.read_stats.num_raw_reads = find_num_reads([x for x in [self.forward_reads, self.reverse_reads] if x is not None])

        if self.read_stats.num_raw_reads == 0:
            logging.error("No reads found for " + self.id + "\n")
            sys.exit(3)
        logging.info("Number of raw reads = " + str(self.read_stats.num_raw_reads) + "\n")
        return

    def prep_reads_for_assembly(self, trimmed_reads: list):
        logging.info("Preparing quality FASTQ files for assembly... ")
        singletons_cat_fq = str(Path(self.output_dir).joinpath("temp_trim").joinpath("singletons.fastq"))
        for fastq in trimmed_reads:
            if re.search(r'pe.[1-2].fq$', fastq):
                self.pe_trimmed.append(fastq)
            else:
                self.se_trimmed.append(fastq)

        try:
            orphans_handler = open(singletons_cat_fq, 'w')
        except IOError:
            logging.error("Unable to open {} for writing.\n".format(singletons_cat_fq))
            sys.exit(3)

        for unpaired_file in self.se_trimmed:
            with open(unpaired_file) as orphan_fq:
                for line in orphan_fq:
                    orphans_handler.write(line)
        orphans_handler.close()
        self.se_trimmed = [singletons_cat_fq]

        logging.info("done.\n")
        return self.pe_trimmed, self.se_trimmed

    def qc_reads(self, background: str, parity: str, adapters: str, executables: dict, num_threads=2):
        filtered_reads = filter_backbone(self, background, executables, parity, num_threads)
        self.read_stats.num_filtered_reads = find_num_reads(filtered_reads)
        self.read_stats.calc_on_target()
        logging.info("{0} reads removed by filtering background ({1}%).\n".format(self.read_stats.num_on_target,
                                                                                  self.read_stats.percent_on_target()))
        # if self.read_stats.num_filtered_reads < 1600:
        #     # The number of reads remaining is too low for assembly (< 20X for a single fosmid)
        #     logging.warning("Number of reads remaining will provide less than 20X coverage for a single fosmid"
        #                     " - skipping this sample\n")
        #     return [], []
        trimmed_reads = quality_trimming(executables["trimmomatic"], self, filtered_reads, parity, adapters, num_threads)
        self.read_stats.num_reads_assembled = find_num_reads(trimmed_reads)
        self.read_stats.num_reads_trimmed = self.read_stats.num_filtered_reads - self.read_stats.num_reads_assembled

        return self.prep_reads_for_assembly(trimmed_reads)

    def prep_nanopore(self, args, executables, nanopore_reads: str):
        raw_nanopore_fasta = read_fasta(nanopore_reads)
        self.read_stats.num_raw_reads = len(raw_nanopore_fasta.keys())
        if self.error_correction is False:
            # Skip the error correction and filter the nanopore reads that don't align to Illumina reads
            # Set nanopore_reads to the filtered raw reads
            self.nanopore = extract_nanopore_for_sample(args, self, executables, raw_nanopore_fasta)
        else:
            self.nanopore = correct_nanopore(args, executables, self)
        # Align the corrected reads to the trim_sequences.fasta file using LAST
        nanopore_background_alignments = align_nanopore_to_background(args, executables, self)
        # Trim the background sequences and reads shorter than 1000bp after removing background
        nanopore_reads, read_stats = filter_background_nanopore(self, nanopore_background_alignments)
        write_trimmed_reads(self, nanopore_reads)
        return

    def get_test_fastq(self):
        try:
            return self.pe_trimmed[0]
        except IndexError:
            return self.se_trimmed[0]

    def assemble_fosmids(self, executables, num_threads):
        k_min, k_max, = determine_k_values(self.get_test_fastq(), self.assembler)
        min_count, cov = determine_min_count(self.read_stats.num_reads_assembled, self.num_fosmids_estimate, k_max)
        with open(Path(self.output_dir).joinpath(f"estimated_coverage.txt"), "w") as f_cov:
            f_cov.write(f"{cov}\n")
        _asm = assemble_fosmids(self, self.assembler, self.assembly_mode, k_min, k_max, min_count, executables, num_threads=num_threads)
        # clean_intermediates(self)
        self.assembled_fosmids = str(Path(self.output_dir).joinpath("temp_assembly").joinpath("final.contigs.fa"))
        if Path(_asm).absolute() != Path(self.assembled_fosmids).absolute():
            os.system(f"mv {_asm} {self.assembled_fosmids}")
        return


# class Miffed(Sample):
#     def __init__(self, sample_id):
#         Sample.__init__(self, sample_id)
#         self.project = ""
#         self.selector = ""
#         self.vector = ""
#         self.screen = ""
#         self.selection = ""
#         self.seq_submission_date = ""
#         self.glycerol_plate = ""
#         self.seq_center = ""
#         self.seq_type = ""
#         self.read_length = ""
#         self.instrument = ""

    # def populate_info(self, fields):
    #     self.project = fields[1]
    #     self.selector = fields[2]
    #     self.vector = fields[3]
    #     self.screen = fields[4]
    #     self.selection = fields[5]
    #     self.seq_submission_date = fields[7]
    #     self.glycerol_plate = fields[8]
    #     self.seq_center = fields[9]
    #     self.seq_type = fields[10]
    #     self.instrument = fields[12]
    #     try:
    #         self.num_fosmids_estimate = int(fields[6])
    #     except ValueError:
    #         logging.error("Number of fosmids field (column 7) is not an integer!\n")
    #         sys.exit(3)
    #     try:
    #         self.read_length = int(fields[11])
    #     except ValueError:
    #         logging.error("Read length field (column 12) is not an integer!\n")
    #         sys.exit(3)

    # def ensure_completeness(self, miffed):
    #     if self.id is None:
    #         logging.warning("Sample name missing in " + miffed + ". Exiting!\n")
    #         sys.exit()
    #     if self.project is None:
    #         logging.warning(self.id + " is not associated with a project!. Skipping this library!\n")
    #         return False
    #     if self.selector is None:
    #         logging.warning(self.id + " is not associated with a human (selector). Skipping this library!\n")
    #         return False
    #     if self.vector is None:
    #         logging.warning(self.id + " is not associated with a vector (e.g. PCC1). Skipping this library!\n")
    #         return False
    #     if self.screen != "in silico" and self.screen != "functional":
    #         logging.warning(self.screen + " was provided for screen in sample " + self.id + ". Is this correct?\n")
    #     return True


class Alignment:
    def __init__(self, fields):
        self.sseqid = fields[0]
        self.slen = int(fields[1])
        self.qseqid = fields[2]
        self.pident = float(fields[3])
        self.length = int(fields[4])
        self.sstrand = fields[5]
        self.sstart = int(fields[6])
        self.send = int(fields[7])
        self.bitscore = float(fields[8])
        self.qstart = int(fields[9])
        self.qend = int(fields[10])
        self.name = ""
        self.direction = ""

    def get_query(self):
        return self.qseqid

    def get_info(self):
        info_string = "contig = " + str(self.sseqid) + \
                      "\tcontig-length = " + str(self.slen) + \
                      "\tclone = " + str(self.name) + \
                      "\tstrand = " + str(self.sstrand) + \
                      "\tdirection = " + str(self.direction) + \
                      "\tstart = " + str(self.sstart) + \
                      "\tend = " + str(self.send) + \
                      "\talignment-length = " + str(self.length) + \
                      "\tpercent-identity = " + str(self.pident)
        return info_string


class EndAlignment(Alignment):
    def parse_fosmid_end_name(self):
        self.direction = "F" if ENDS_FW_FLAG in self.qseqid else "R"
        self.name = get_fosmid_end_name(self.qseqid)


class MyFormatter(logging.Formatter):

    error_fmt = "%(levelname)s - %(module)s, line %(lineno)d:\n%(message)s"
    warning_fmt = "%(levelname)s:\n%(message)s"
    debug_fmt = "%(asctime)s\n%(message)s"
    info_fmt = "%(message)s"

    def __init__(self):
        super().__init__(fmt="%(levelname)s: %(message)s",
                         datefmt="%d/%m %H:%M:%S")

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._style._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._style._fmt = MyFormatter.debug_fmt

        elif record.levelno == logging.INFO:
            self._style._fmt = MyFormatter.info_fmt

        elif record.levelno == logging.ERROR:
            self._style._fmt = MyFormatter.error_fmt

        elif record.levelno == logging.WARNING:
            self._style._fmt = MyFormatter.warning_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._style._fmt = format_orig

        return result


def subprocess_helper(cmd_list, collect_all=True, graceful=False):
    """
    Wrapper function for opening subprocesses through subprocess.Popen()

    :param cmd_list: A list of strings forming a complete command call
    :param collect_all: A flag determining whether stdout and stderr are returned
    via stdout or just stderr is returned leaving stdout to be written to the screen
    :param graceful: Allows an executable to exit with a non-zero returncode and still return
    :return: A string with stdout and/or stderr text and the returncode of the executable
    """
    stdout = ""
    if collect_all:
        proc = subprocess.Popen(' '.join(cmd_list),
                                shell=True,
                                preexec_fn=os.setsid,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
        stdout = proc.communicate()[0].decode("utf-8")
    else:
        proc = subprocess.Popen(' '.join(cmd_list),
                                shell=True,
                                preexec_fn=os.setsid)
        proc.wait()

    # Ensure the command completed successfully
    if proc.returncode != 0 and not graceful:
        logging.error(cmd_list[0] + " did not complete successfully! Command used:\n" +
                      ' '.join(cmd_list) + "\nOutput:\n" + stdout)
        sys.exit(19)

    return stdout, proc.returncode


def bam2fastq(input_file: str, prefix: str, output_dir: str, samtools_exe: str, parity="pe"):
    """
    :param input_file: A BAM file with the aligned reads
    :param sample_id: Name of the sample - to be used for identifying output reads
    :param output_dir: Path to write the reads
    :param samtools_exe: Path to the samtools executable, version
    :param parity: Argument indicating whether the reads are from a paired-end (pe) or single-end (se) library
    :return: List of paths to FASTQ files converted from the BAM input_file
    """
    output_dir = os.path.abspath(output_dir)
    fastq_extract = [samtools_exe, "fastq", "--verbosity", str(1), "-N"]
    if parity == "pe":
        fastq_files = [os.path.join(output_dir, prefix + ".1.fastq"),
                       os.path.join(output_dir, prefix + ".2.fastq")]
        fastq_extract += ["-1", fastq_files[0]]
        fastq_extract += ["-2", fastq_files[1]]
        # Need to include the singletons file so orphaned reads are removed
        fastq_extract += ["-s", os.path.join(output_dir, prefix + ".singletons.fastq")]
        fastq_extract += [input_file]
    elif parity == "se":
        fastq_files = [os.path.join(output_dir, prefix + ".fastq")]
        fastq_extract += [input_file, ">", fastq_files[0]]
    else:
        raise AssertionError("Unknown value for parity '{}'. Expected either 'pe' or 'se'.\n")
    
    bam2fastq_stdout = os.path.join(output_dir, "samtools_fastq.stdout")
    logging.info("Converting BAM file to FastQ format... ")
    stdout, retcode = subprocess_helper(fastq_extract, True)
    logging.info("done.\n")

    with open(bam2fastq_stdout, 'w') as stdout_file:
        stdout_file.write(stdout)

    for fq_path in fastq_files:
        if not os.path.isfile(fq_path):
            logging.error("Path to 'samtools fastq' output file {} doesn't exist.\n"
                          "Command used:\n"
                          "{}".format(fq_path, ' '.join(fastq_extract)))
            sys.exit(5)

    return fastq_files

ENDS_NAME_REGEX = None
ENDS_FW_FLAG = None
def get_options(sys_args: list):
    parser = argparse.ArgumentParser(
        prog="fabfos",
        description="Pipeline for filtering, assembling and organizing fosmid sequence information.\n",
        add_help=False
    )
    reqs = parser.add_argument_group(title="Required arguments")
    reads = parser.add_argument_group(title="Sequence read-specific arguments")
    nanopore = parser.add_argument_group(title="Nanopore-specific [development] options")
    opts = parser.add_argument_group(title="Optional arguments")
    misc = parser.add_argument_group(title="Miscellaneous options")

    reqs.add_argument("-b", "--background", type=str, required=False,
                      help="Path to fosmid background fasta")
    reqs.add_argument("--vector", type=str, required=False,
                      help="Path to vector backbone fasta, used for estimating pool size, required if --pool-size not given")
    reqs.add_argument("-n", "--pool-size", type=int, required=False,
                      help="Estimate of number of fosmids in pool, required if --vector not given")
    
    nanopore.add_argument("--nanopore_reads", type=str, default=None, required=False,
                          help="Experimental feature! A FASTA file containing nanopore reads to be used in assembly.")
    nanopore.add_argument("--skip_correction", action="store_true", default=False, required=False,
                          help="Experimental feature! Do not perform error-correction of nanopore reads using proovread")

    reads.add_argument("-r", "--reads", type=str, required=False,
                       help="Path to the forward strand file or the interleaved paired-end file. Can be in either FastQ or BAM format.")
    reads.add_argument("-2", "--reverse", type=str, required=False,
                       help="Path to the reverse-end read file (if applicable)")
    reads.add_argument("-t", "--type", choices=["B", "F"], required=False, default="F",
                       help="Enter B if input type is BAM, F for FastQ. [ DEFAULT = 'F' ]")
    reads.add_argument("-i", "--interleaved", required=False, default=False, action="store_true",
                       help="Flag indicating the reads are interleaved "
                            "(i.e. forward and reverse pairs are in the same file). [DEFAULT = False]")
    reads.add_argument("-p", "--parity", required=False, default="pe", choices=["pe", "se"],
                       help="Specifying the sequencing chemistry used, either paired-end (pe) or single-end (se). "
                            "[DEFAULT = 'pe']")

    opts.add_argument("-a", "--assembler", required=False,
                      choices=["spades_meta", "spades_isolate", "spades_sc", "megahit"], default="megahit",
                      help="Genome assembly software to use. [DEFAULT = megahit]")
    opts.add_argument("-o", "--output", type=str, required=False,
                      default=os.path.abspath("./"),
                      help="path to temp. workspace [DEFAULT = /tmp]")
    opts.add_argument("-T", "--threads", type=str, required=False, default=str(8),
                      help="The number of threads that can be used [DEFAULT = 8]")
    opts.add_argument("-e", "--ends", required=False, default=None,
                      help="FASTA file containing fosmid ends - these will be used for alignment to fosmid contigs.")
    opts.add_argument("--ends-name-regex", required=False, default=None,
                      help='regex for getting name of endseq., ex. "\\w+_\\d+" would get ABC_123 from ABC_123_FW')
    opts.add_argument("--overwrite", required=False, default=False, action="store_true",
                      help="overwrite the output directory if it already exists")
    opts.add_argument("--ends-fw-flag", required=False, default=None,
                      help='string that marks forward ends, ex. "_FW" would cause ABC_123_FW and ABC_123_RE to be assigned to forward and reverse respectively.')
    opts.add_argument("-l", "--min-length", type=int, required=False, default=1000,
                      help="Assembled contigs shorter than this will be discarded [DEFAULT = 1000]")
    
    misc.add_argument("-v", "--version", default=False, action="store_true",
                      help="Print the FabFos version and exit.")
    misc.add_argument("--verbose", required=False, default=False, action="store_true",
                      help="Increase the level of verbosity in runtime log.")
    misc.add_argument("-h", "--help",
                      action="help", help="Show this help message and exit")

    args = parser.parse_args(sys_args)

    if args.ends:
        # yes, using global variables is lazy, but this is also a 2k+ line python file...
        name_pattern = args.ends_name_regex
        assert name_pattern is not None, f"if using endseq, also provide name pattern"
        fw_flag = args.ends_fw_flag
        assert fw_flag is not None, f"if using endseq, also provide forward direction flag"
        global ENDS_NAME_REGEX, ENDS_FW_FLAG
        ENDS_NAME_REGEX = name_pattern
        ENDS_FW_FLAG = fw_flag
    return args


def os_type():
    """
    Return the operating system of the user
    """
    x = sys.platform
    if x:

        hits = re.search(r'darwin', x, re.I)
        if hits:
            return 'mac'

        hits = re.search(r'win', x, re.I)
        if hits:
            return 'win'

        hits = re.search(r'linux', x, re.I)
        if hits:
            return 'linux'


def prep_logging(log_file_name: str, verbose=False):
    logging.basicConfig(level=logging.DEBUG,
                        filename=log_file_name,
                        filemode='w',
                        datefmt="%d/%m %H:%M:%S",
                        format="%(asctime)s %(levelname)s:\n%(message)s")
    if verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.INFO

    # Set the console handler normally writing to stdout/stderr
    ch = logging.StreamHandler()
    ch.setLevel(logging_level)
    ch.terminator = ''

    formatter = MyFormatter()
    ch.setFormatter(formatter)
    logging.getLogger('').addHandler(ch)
    return

def review_arguments(args, fabfos: FabFos) -> None:
    """
    Function that ensures the proper information has been provided

    :param args: parsed command-line arguments from get_options()
    :param fabfos: Instance of the FabFos class
    :return: None
    """

    # Review the provided arguments:
    if not args.reads:
        logging.error("the following argument is required: -r/--reads\n")
        sys.exit(1)
    if not os.path.isfile(args.reads):
        logging.error(args.reads + " is not a file\n")
        sys.exit(3)
    if args.reverse and not os.path.isfile(args.reverse):
        logging.error(args.reverse + " is not a file\n")
        sys.exit(3)

    # warning: nanopore is experimental
    if args.nanopore_reads is not None:
        logging.warn(f"using nanopore reads is an experimental feature!")

    # at least one of --vector or --pool-size must be provided
    if args.vector is None and args.pool_size is None:
        logging.error("at least one of --vector or --size must be provided\n")
        sys.exit(6)

    # Check if the background (FASTA including vector and source genome) file was provided and exists
    if not args.background:
        logging.error("a fasta for the host background is required: -b/--background\n")
        sys.exit(1)
    if not os.path.isfile(args.background):
        logging.error(args.background + " does not exist!\n")
        sys.exit(3)

    if args.reverse and args.interleaved:
        logging.error("Reads cannot be interleaved and also have separate forward- and reverse-FASTQ files!\n")
        sys.exit(3)

    if args.parity == "se" and args.type == "B":
        logging.error("FabFos is not capable of dealing with BAM files generated from single-end sequencing.\n"
                      "Please convert the BAM to a FASTQ file, using samtools or bam2fastq, and provide those.\n")
        sys.exit(5)
    return


def is_exe(program):
    return os.path.isfile(program) and os.access(program, os.X_OK)


def which(program):
    f_path, f_name = os.path.split(program)
    if f_path:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            prefix, ext = os.path.splitext(exe_file)
            if is_exe(exe_file):
                return exe_file
            elif ext == ".jar" and os.path.isfile(exe_file):
                return exe_file
    return None


def find_executables(assembler: str, nanopore=False) -> dict:
    """
    Find the executables that are to be used in the pipeline:
    BWA, MEGAHIT, Trimmomatic, samtools

    :return: Dictionary with executable paths indexed by executable names
    """
    required_execs = ["bwa", "samtools", "trimmomatic", "blastn", "makeblastdb"]
    np_execs = ["proovread", "canu"]
    executables = dict()
    for executable in required_execs:
        executables[executable] = which(executable)
        if executables[executable] is None:
            raise EnvironmentError("Unable to find executable for " + executable)

    if assembler.startswith("spades"):
        executables["spades"] = which("spades.py")
    elif assembler == "megahit":
        executables["megahit"] = which("megahit")
    else:
        logging.error("Unsupported assembler specified: '{}'.\n".format(assembler))
        sys.exit(3)

    if nanopore:
        for executable in np_execs:
            executables[executable] = which(executable)
            if executables[executable] is None:
                raise EnvironmentError("Unable to find executable for " + executable)
    return executables


def find_dependency_versions(exe_dict: dict) -> dict:
    """
    Function for retrieving the version numbers for each executable in exe_dict

    :param exe_dict: A dictionary mapping names of software to the path to their executable
    :return: A formatted string with the executable name and its respective version found
    """
    versions_dict = dict()

    simple_v = ["megahit", "spades"]
    version_param = ["blastn", "makeblastdb"]
    no_params = ["bwa", "samtools", ]
    for exe in exe_dict:
        ##
        # Get the version statement for the software
        ##
        versions_dict[exe] = ""
        if exe in simple_v:
            stdout, returncode = subprocess_helper([exe_dict[exe], "-v"])
        elif exe in version_param:
            stdout, returncode = subprocess_helper([exe_dict[exe], "-version"])
        elif exe in no_params:
            stdout, returncode = subprocess_helper([exe_dict[exe]], graceful=True)
        elif exe == "trimmomatic":
            stdout, returncode = subprocess_helper([exe_dict[exe], "-version"])
            versions_dict[exe] = stdout.strip()
            continue
        else:
            logging.warning("Unknown version command for " + exe + ".\n")
            continue
        ##
        # Identify the line with the version number (since often more than a single line is returned)
        ##
        for line in stdout.split("\n"):
            matches = re.findall(r"\d+\.\d+(?:\.[^\s,\.]+)?", line)
            if len(matches) == 0: continue
            versions_dict[exe] = matches[0]
        if not versions_dict[exe]:
            logging.debug("Unable to find version for " + exe + ".\n")

    return versions_dict


def validate_dependency_versions(dep_versions: dict) -> bool:
    """
    Ensure the versions of executables are compatible with FabFos.

    :param dep_versions: A dictionary mapping executable names to the installed version
    :return: Boolean indicating whether all of the installed executable dependencies (not Python libraries) are the
    right version.
    """
    reqd_versions = {"samtools": "1.10",
                     "trimmomatic": "0.39"}
    for dep, min_v in reqd_versions.items():
        if version.parse(dep_versions[dep]) < version.parse(min_v):
            logging.warning("{} version found ('{}') is not compatible with FabFos - {} or later required.\n"
                            "".format(dep, dep_versions[dep], reqd_versions[dep]))
            return False

    return True


def summarize_dependency_versions(dep_versions: dict) -> None:
    versions_string = "FabFos v{}".format(VERSION)
                    #   "Software versions used:\n"
    # ##
    # # Format the string with the versions of all software
    # ##
    # for exe in sorted(dep_versions):
    #     n_spaces = 20 - len(exe)
    #     versions_string += "\t" + exe + ' ' * n_spaces + dep_versions[exe] + "\n"

    if logging.getLogger().hasHandlers():
        logging.info(versions_string)
    else:
        print(versions_string)
    return


def check_index(path, bwa_path):
    extensions = ['.bwt', '.pac', '.ann', '.amb', '.sa']
    index_parts = [path+e for e in extensions]
    for part in index_parts:
        if not os.path.isfile(part):
            # logging.info("\nUnable to find BWA index for " + path + ". Indexing now... ")
            logging.info("\nIndexing "+Path(path).name+"...")
            index_command = [bwa_path, "index", path]
            index_dir = os.path.dirname(os.path.abspath(path))
            index_command += ["1>", index_dir + os.sep + "bwa_index.stdout"]
            index_command += ["2>", index_dir + os.sep + "bwa_index.stderr"]

            subprocess_helper(index_command, graceful=False)

            logging.info("done.\nResuming alignment... ")
            break
    return


def bwa_mem_wrapper(bwa_exe, index, fwd_fq, output_dir,
                    rev_fq=None, num_threads=2, interleaved=False) -> str:
    align_command = [bwa_exe, "mem",
                     "-t", str(num_threads),
                     index,
                     fwd_fq]
    if rev_fq is not None:
        align_command.append(rev_fq)
    if interleaved:
        align_command.append("-p")

    align_command.append("1>")
    sam_file = str(Path(output_dir).joinpath("aligned.sam"))
    align_command.append(sam_file)
    align_command += ["2>", "/dev/null"]

    subprocess_helper(cmd_list=align_command, graceful=False)

    return sam_file


def filter_backbone(sample: Sample, background: str, executables: dict, parity: str, num_threads=2):
    """
    Function to generate fastq files that do not contain any sequences in `background`.
    Depends on bwa, samtools, and bam2fastq

    :param sample: Miffed object with information of current sample
    :param background: Path to a FASTA file containing the genomic sequences to be removed
    :param executables: Dictionary containing paths to executables that is indexed by the executable name
    :param parity: Argument indicating whether the reads are from a paired-end (pe) or single-end (se) library
    :param num_threads: Number of threads available for BWA and samtools to use
    :return: list of fastq files containing the filtered reads
    """
    logging.info("Filtering off-target reads... ")
    folder = str(Path(sample.output_dir).joinpath("temp_filter"))
    os.makedirs(folder, exist_ok=True)
    shutil.copy(background, folder)
    bpath = f"{Path(folder).joinpath(Path(background).name)}"
    check_index(bpath, executables["bwa"])
    sam_file = bwa_mem_wrapper(bwa_exe=executables["bwa"], index=bpath, num_threads=num_threads,
                               fwd_fq=sample.forward_reads, rev_fq=sample.reverse_reads,
                               output_dir=folder, interleaved=sample.interleaved)

    # Use samtools to convert sam to bam
    bam_convert = [executables["samtools"], "view", "-bS", "-f", str(4), "-@", str(num_threads), sam_file, "|"]
    # Piping to samtools sort
    bam_file = folder + os.sep  + "sorted.bam"
    bam_convert.append(executables["samtools"])
    bam_convert += ["sort", "-@", str(num_threads), "-"]
    bam_convert += ["-o", bam_file]
    p_samtools_stdout = open(folder + os.sep + "samtools.stdout", 'w')

    stdout, retcode = subprocess_helper(cmd_list=bam_convert, graceful=False)
    p_samtools_stdout.write(stdout)
    logging.info("done.\n")
    p_samtools_stdout.close()

    # extract the unaligned reads using bam2fastq
    filtered_reads = bam2fastq(bam_file, "filtered", folder, executables["samtools"], parity)
    return filtered_reads


def get_reference_names_from_sam(sam_file):
    reference_names = set()

    try:
        sam = open(sam_file, 'r')
    except IOError:
        logging.error("Unable to open " + sam_file + " for reading!\n")
        sys.exit(3)

    line = sam.readline()
    # Skip the header lines:
    while line[0] == "@":
        line = sam.readline()
    # Read the alignments:
    while line:
        fields = line.strip().split("\t")
        if fields[2] != "*":
            reference_names.add(fields[2])
        line = sam.readline()

    sam.close()
    return reference_names

def write_new_fasta(fasta_dict: dict, fasta_name: str, headers=None):
    """
    Function for writing sequences stored in dictionary to file in FASTA format; optional filtering with headers list

    :param fasta_dict: A dictionary containing headers as keys and sequences as values
    :param fasta_name: Name of the FASTA file to write to
    :param headers: Optional list of sequence headers. Only fasta_dict keys in headers will be written
    :return:
    """
    headers = list(headers)
    try:
        fa_out = open(fasta_name, 'w')
    except IOError:
        raise IOError("Unable to open " + fasta_name + " for writing!")

    for name in fasta_dict.keys():
        seq = fasta_dict[name]
        if headers is None:
            fa_out.write(name + "\n")
            fa_out.write(seq + "\n")
        elif name[1:] in headers:
            fa_out.write(name + "\n")
            fa_out.write(seq + "\n")

    fa_out.close()
    return


def extract_nanopore_for_sample(args, sample: Sample, executables: dict, raw_nanopore_fasta: dict):
    """
    Function aligns the Illumina reads to long reads and creates a new FASTA file of reads that were aligned to

    :param args:
    :param sample:
    :param executables:
    :param raw_nanopore_fasta:
    :return: Name of the output FASTA file
    """
    logging.info("Aligning Illumina reads to long reads... ")

    folder = str(Path(sample.output_dir).joinpath("temp_nanopore"))
    select_nanopore = folder + os.sep + "select_nanopore.fasta"
    check_index(args.nanopore_reads, executables["bwa"])

    filtered_forward = ""
    filtered_reverse = ""
    for fastq in sample.pe_trimmed:
        if re.search(r'pe.1.fq$', fastq):
            filtered_forward = fastq
        else:
            filtered_reverse = fastq

    align_command = [executables["bwa"], "mem", "-t", args.threads, args.nanopore_reads, filtered_forward]
    if args.reverse:
        align_command.append(filtered_reverse)
    align_command.append("1>")
    sam_file = folder + os.sep + "nanopore_hits.sam"
    align_command.append(sam_file)
    align_command += ["2>", "/dev/null"]

    subprocess_helper(align_command, False)

    logging.info("done.\n")

    # Parse the SAM file and find the names of the nanopore reads that were mapped to
    logging.info("Parsing SAM file for names of the long reads that were mapped to... ")
    mapped_long_reads = get_reference_names_from_sam(sam_file)
    os.remove(sam_file)
    logging.info("done.\n")

    logging.info("Writing sequences to " + select_nanopore + "... ")
    write_new_fasta(raw_nanopore_fasta, select_nanopore, mapped_long_reads)
    logging.info("done.\n")

    return select_nanopore


def correct_nanopore(args, executables, sample):
    """
    A function to correct errors in the nanoopore reads using short, accurate Illumina reads and proovread

    :param sample: Miffed object with information of current sample
    :param args: command-line arguments list
    :param executables: A dictionary containing names of executables as keys and their respective file paths as values
    :return: Name of the corrected and trimmed nanopore reads
    """
    logging.info("Correcting errors in " + args.nanopore_reads + " using proovread... ")

    folder = str(Path(sample.output_dir).joinpath("temp_nanopore"))

    proovread_prefix = folder + os.sep + f"{sample.id}_proovread."
    proovread_command = [executables["proovread"], "--threads", args.threads]
    proovread_command += ["--long-reads=" + args.nanopore_reads]
    proovread_command += ["--short-reads=" + sample.forward_reads]
    if args.reverse:
        proovread_command += ["--short-reads=" + sample.reverse_reads]
    proovread_command += ["-p", folder]  # prefix to output files
    try:
        correct_stdout = open(proovread_prefix + "stdout", 'w')
    except IOError:
        raise IOError("ERROR: cannot open file: " + proovread_prefix + "stdout")

    stdout, retcode = subprocess_helper(cmd_list=proovread_command)

    correct_stdout.write(stdout)
    correct_stdout.close()
    logging.info("done.\n")
    return folder + os.sep + "proovread.trimmed.fa"


def quality_trimming(trimmomatic_exe: str, sample: Sample, filtered_reads: list, parity: str, adapters: str,
                     num_threads=2) -> list:
    """
    Wrapper for trimmomatic

    :param sample: Miffed object with information of current sample
    :param filtered_reads: list of background-filtered fastq files
    :param parity: String indicating the sequencing chemistry library for the library [pe (default)|se]
    :param adapters: Path to a file containing Illumina adapter sequences
    :param trimmomatic_exe: Paths to the Trimmomatic executable
    :param num_threads: Number of threads available parallel processing
    :return: list of quality-trimmed fastq files
    """
    logging.info("Trimming reads... ")

    folder = str(Path(sample.output_dir).joinpath("temp_trim"))
    os.makedirs(folder, exist_ok=True)
    trim_prefix = folder + os.sep

    trimmomatic_command = [trimmomatic_exe]
    if parity == "pe":
        trimmomatic_command.append("PE")
        trimmomatic_outputs = [trim_prefix + "pe.1.fq", trim_prefix + "se.1.fq",
                               trim_prefix + "pe.2.fq", trim_prefix + "se.2.fq"]
    elif parity == "se":
        trimmomatic_command.append("SE")
        trimmomatic_outputs = [trim_prefix + "se.1.fq"]
    else:
        logging.error("Unexpected read-level parity: '{}'. Either 'pe' or 'se' are expected.\n".format(parity))
        sys.exit(17)
    trimmomatic_command += ["-threads", num_threads]
    trimmomatic_command += filtered_reads
    trimmomatic_command += trimmomatic_outputs

    adapters = "ILLUMINACLIP:" + adapters + "TruSeq3-PE.fa:2:3:10"
    trimmomatic_command.append(adapters)
    trimmomatic_command += ["LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"]
    try:
        trim_stdout = open(trim_prefix + "stdout.txt", 'w')
    except IOError:
        logging.error("Cannot open file: " + trim_prefix + "stdout.txt for writing.\n")
        sys.exit(3)

    stdout, retcode = subprocess_helper(cmd_list=trimmomatic_command)

    trim_stdout.write(stdout)
    trim_stdout.close()
    logging.info("done.\n")
    return trimmomatic_outputs


def determine_k_values(test_fastq: str, assembler: str):
    k_min = 71
    if assembler == "spades":
        k_max = 127
    else:
        k_max = 241
    max_read_length = 0
    read_lengths = list()
    # Sample the first paired-end FASTQ file

    x = 0
    for _, seq in read_fastq(test_fastq):
        read_length = len(seq)
        read_lengths.append(read_length)
        if read_length > max_read_length:
            max_read_length = read_length
        if x > 1E5:
            break
        x += 1

    if max_read_length < k_max:
        if max_read_length % 2 == 0:
            k_max = max_read_length - 3
        else:
            k_max = max_read_length - 2
        if max_read_length < k_min:
            k_min = k_max
    return k_min, k_max


def spades_wrapper(sample: Sample, k_min: int, k_max: int, min_count: int,
                   spades_exe=None, mode="isolate", num_threads=2) -> (str, str, str):
    # The following is for assembling with SPAdes:
    if not spades_exe:
        spades_exe = "spades.py"

    spades_command = [spades_exe]
    if sample.parity == "pe":
        fwd_fq, rev_fq = sample.retrieve_separate_pe_fq()
        spades_command += ["-1", fwd_fq]
        if rev_fq:
            spades_command += ["-2", rev_fq]
    if len(sample.se_trimmed) >= 1:
        spades_command += ["-s", ','.join(sample.se_trimmed)]
    spades_command += ["--memory", str(20)]
    spades_command += ["--threads", str(num_threads)]
    spades_command += ["-k", ','.join([str(k) for k in range(k_min, k_max, 10)])]
    spades_command += ["-o", sample.output_dir + "temp_assembly"]
    spades_command.append("--{}".format(mode))
    if mode != "meta":
        spades_command += ["--cov-cutoff", str(min_count)]

    subprocess_helper(spades_command)

    contigs = sample.output_dir + "temp_assembly" + os.sep + "contigs.fasta"
    asm_log = sample.output_dir + "temp_assembly" + os.sep + "spades.log"
    scaffolds = sample.output_dir + "temp_assembly" + os.sep + "scaffolds.fasta"
    return contigs, scaffolds, asm_log


def megahit_wrapper(sample: Sample, megahit_exe: str, k_min: int, k_max: int, min_count: int, num_threads=2):
    megahit_command = [megahit_exe]
    if sample.parity == "pe":
        fwd_fq, rev_fq = sample.retrieve_separate_pe_fq()
        megahit_command += ["-1", fwd_fq]
        if rev_fq:
            megahit_command += ["-2", rev_fq]
    if len(sample.se_trimmed) >= 1:
        megahit_command += ["--read", ','.join(sample.se_trimmed)]
    megahit_command += ["--k-min", str(k_min), "--k-max", str(k_max)]
    megahit_command += ["--min-count", str(min_count)]
    megahit_command += ["--k-step", str(10)]
    megahit_command += ["--memory", str(0.25)]
    megahit_output_dir = sample.output_dir + "temp_assembly"
    megahit_command += ["--out-dir", megahit_output_dir]
    megahit_command += ["--num-cpu-threads", str(num_threads)]
    megahit_command.append("--no-mercy")  # Recommended for high-coverage datasets

    subprocess_helper(megahit_command)

    contigs = sample.output_dir + "temp_assembly" + os.sep + "final.contigs.fa"
    asm_log = sample.output_dir + "temp_assembly" + os.sep + "log"
    scaffolds = ""
    return contigs, scaffolds, asm_log


def assemble_fosmids(sample: Sample, assembler: str, assembly_mode: str, k_min: int, k_max: int, min_count: int,
                     executables: dict, k_step=10, num_threads=2):
    """
    Wrapper function for the assembly process - multi-sized de Bruijn graph based assembler for metagenomes

    :param sample: Miffed object with information of current sample
    :param assembler: String indicating which assembler [megahit|spades] should be used
    :param assembly_mode: When assembling with SPAdes, the assembly mode can be toggled between 'isolate' or 'meta'
    :param k_min: 71; painstakingly determined to be optimal for fosmids sequenced with Illumina
    :param k_max: 241; painstakingly determined to be optimal for fosmids sequenced with Illumina
    :param min_count: The minimum k-mer abundance to be used by megahit for building the succinct DBG
    :param executables: Dictionary containing paths to executables that is indexed by the executable name
    :param num_threads: Number of threads available for parallel computation
    :return:
    """
    asm_param_stmt = "Assembling sequences using {}\n".format(assembler)
    # if assembly_mode != "NA":
    #     asm_param_stmt += "Assembly mode is '{}'\n".format(assembly_mode)
    # asm_param_stmt += "Parameters:\n--k-min = {}\t--k-max = {}\t--k-step = {}".format(k_min, k_max, k_step)
    # if assembly_mode != "meta":
    #     asm_param_stmt += "\t--min-count = {}\n".format(min_count)
    # else:
    #     asm_param_stmt += "\n"
    logging.info(asm_param_stmt)

    if assembler == "spades":
        f_contigs, f_scaffolds, asm_log = spades_wrapper(sample=sample, mode=assembly_mode,
                                                         k_min=k_min, k_max=k_max, min_count=min_count,
                                                         num_threads=num_threads)
    elif assembler == "megahit":
        f_contigs, f_scaffolds, asm_log = megahit_wrapper(sample=sample, megahit_exe=executables["megahit"],
                                                          k_min=k_min, k_max=k_max, min_count=min_count,
                                                          num_threads=num_threads)
    else:
        logging.error("Unknown assembly software '{}' requested.\n".format(assembler))
        sys.exit(3)

    if not os.path.isfile(f_contigs):
        logging.error("{0} assembly outputs were not created!\n"
                      "Check {1} log for an error.".format(f_contigs, sample.output_dir + "assembly" + os.sep))
        sys.exit(3)

    # logging.info("Cleaning up assembly outputs... ")
    # Move the output files to their final destinations
    # shutil.move(f_contigs, sample.output_dir + os.sep + sample.id + "_contigs.fasta")
    # shutil.move(asm_log, sample.output_dir + os.sep + assembler + "_log.txt")
    # if len(f_scaffolds) > 0:
    #     # Only available for SPAdes
    #     shutil.move(f_scaffolds,
    #                 sample.output_dir + os.sep + sample.id + "_scaffolds.fasta")
    # shutil.rmtree(sample.output_dir + "assembly" + os.sep)
    logging.info("Assembly done.\n")
    return f_contigs


def read_fastq_to_dict(fastq_file: str) -> dict:
    fq_dict = {}
    acc = 0
    matepair_re = re.compile(r".*/[12]$")
    logging.info("Reading FASTQ file '{}'... ")
    for name, seq in read_fastq(fastq_file):
        if not matepair_re.match(name):
            if acc % 2:
                name += "/2"
            else:
                name += "/1"
        fq_dict[name] = seq
        acc += 1

    logging.info("done.\n")
    return fq_dict


def deinterleave_fastq(fastq_file: str, output_dir) -> (str, str):
    # TODO: support gzipped files
    logging.info("De-interleaving forward and reverse reads in " + fastq_file + "... ")

    if not output_dir:
        output_dir = os.path.dirname(fastq_file)
    prefix = '.'.join(os.path.basename(fastq_file).split('.')[:-1])
    fwd_fq = output_dir + os.sep + prefix + "_R1.fastq"
    rev_fq = output_dir + os.sep + prefix + "_R2.fastq"

    # https://gist.github.com/nathanhaigh/3521724
    os.system(f"""\
        {HERE}/deinterleave_fastq.sh <{fastq_file} {fwd_fq} {rev_fq}
    """)

    logging.info("done.\n")
    return fwd_fq, rev_fq


def match_file_with_sample(sample_name, file_list):
    regex_sample = re.compile(r'^' + re.escape(sample_name) + r"[._-].*")
    matches = []
    for f_name in file_list:
        if regex_sample.match(os.path.basename(f_name)):
            matches.append(os.path.join(os.getcwd(), f_name))
        else:
            # File is not of the current sample
            pass
    # Ensure there is a single fastq file for this sample
    if len(matches) > 1:
        logging.error("More than two files match sample ID '{}' in reads directory: {}\n"
                      "Please concatenate files from the same library and re-run.\n".format(sample_name,
                                                                                            ', '.join(matches)))
        sys.exit(17)
    elif len(matches) == 0:
        logging.error("Unable to find file matching sample ID '{}'".format(sample_name))
        sys.exit(19)

    return matches.pop()


def find_raw_reads(reads_dir: str, sample_id: str, parity="pe", reverse="") -> dict:
    """
    finds the raw read files in the path and stores the file names in a dictionary.

    :param reads_dir: Path to a directory containing interleaved or forward-end FASTQ files
    :param sample_id: Name of the sample, used for finding the corresponding FASTQ files
    :param parity: String indicating the sequencing chemistry library for the library [pe (default)|se]
    :param reverse: Path to a directory containing reverse-oriented FASTQ files (optional)
    :return: A dictionary divided by forward and reverse sense
    """
    # Prep the dictionary for storing matched files
    raw_reads = dict()
    raw_reads["forward"] = ""
    raw_reads["reverse"] = ""
    forward_reads = []
    reverse_reads = []

    # The list of extensions supported by FabFos
    extensions = ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]
    wildcards = [os.path.join(reads_dir, ext) for ext in extensions]
    for file_type in wildcards:
        forward_reads.extend(glob.glob(file_type))

    if len(forward_reads) == 0:
        logging.error("Unable to locate fastq files. File extensions must be either 'fastq' or 'fq'\n")
        sys.exit(3)
    raw_reads["forward"] = match_file_with_sample(sample_id, forward_reads)
    if reverse and parity == "pe":
        wildcards = [os.path.join(reverse, ext) for ext in extensions]
        for file_type in wildcards:
            reverse_reads.extend(glob.glob(file_type))
        raw_reads["reverse"] = match_file_with_sample(sample_id, reverse_reads)

    # Check to make sure there are fastq files
    elif len(raw_reads.values()) == 0:
        logging.error("Unable to locate fastq files " + sample_id + "\n")
    return raw_reads


def get_exclude_input(sample_id):
    response = input("Should this sample be excluded from analysis? [y|n]")
    if response == "y":
        exclude = True
        logging.info("Excluding " + str(sample_id) + "\n")
    elif response == "n":
        exclude = False
    else:
        logging.info("Unknown response" + "\n")
        exclude = get_exclude_input(sample_id)
    return exclude


def get_overwrite_input(sample_id):
    response = input("Output directory for " + sample_id + " already exists with files. Overwrite? [y|n]")
    if response == "y":
        overwrite = True
        logging.info("Overwriting " + str(sample_id) + "\n")
    elif response == "n":
        overwrite = False
    else:
        logging.info("Unknown response" + "\n")
        overwrite = get_overwrite_input(sample_id)
    return overwrite


# def get_assemble_input(sample):
#     response = input("Assembled fosmids for " + sample.id + " already exists. Do you want to re-assemble? [y|n]")
#     if response == "y":
#         assemble = True
#         logging.info("Re-assembling " + str(sample.id) + "\n")
#     elif response == "n":
#         assemble = False
#         logging.info("Using " + str(sample.assembled_fosmids) + "\n")
#     else:
#         logging.info("Unknown response" + "\n")
#         assemble = get_overwrite_input(sample.id)
#     return assemble

# due to horrific design, this function attempts to sync the states between the passed in objects
def ghetto_sync_args(args, fabfos: FabFos, sample: Sample):
    sample.output_dir = os.sep.join([fabfos.workspace]) + os.sep
    if args.nanopore_reads:
        sample.nanopore = True
    if args.skip_correction is True:
        sample.error_correction = False
    sample.parity = args.parity
    sample.interleaved = args.interleaved
    # Split args.assembler into assembler and mode attributes
    if args.assembler.startswith("spades"):
        sample.assembler, sample.assembly_mode = args.assembler.split('_')
    else:
        sample.assembler, sample.assembly_mode = args.assembler, "NA"

def clean_intermediates(sample):
    """
    Function removes largely useless alignment files

    :param sample: Miffed object with information of current sample
    :return:
    """
    # Remove the alignment files
    working_dir = sample.output_dir + os.sep
    os.remove(working_dir + sample.id + ".sam")
    os.remove(working_dir + sample.id + "_sorted.bam")

    # Remove bam2fastq outputs:
    for std_file in glob.glob(working_dir + "bam2fastq*"):
        os.remove(std_file)

    # Remove trimmomatic outputs
    trimmomatic_outputs = glob.glob(working_dir + sample.id + "_trim*fq")
    for trim_out in trimmomatic_outputs:
        os.remove(trim_out)

    os.remove(working_dir + "singletons.fastq")
    return


def find_num_reads(file_list: list) -> int:
    """
    Function to count the number of reads in all FASTQ files in file_list

    :param file_list: A list of FASTQ files
    :return: integer representing the number of reads in all FASTQ files provided
    """
    logging.info("Finding number of reads in fastq files... ")
    num_reads, nucleotides = 0, 0
    for f in file_list:
        proc = subprocess.Popen(f"cat {f} | awk 'NR % 4 == 2' | wc -cl",
                        shell=True,
                        preexec_fn=os.setsid,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT)
        stdout = proc.communicate()[0].decode("utf-8")
        toks = stdout[:-1].strip()
        if "\t" in toks: toks = toks.split("\t")
        else: toks = [t for t in toks.split(" ") if len(t)>0]
        nr, nuc = [int(t) for t in toks]
        num_reads += nr
        nucleotides += nuc
    logging.info("done.\n")
    return num_reads


def map_ends(executables, ends, sample: Sample):
    """
    Function for using blastn to align fosmid ends in a file to find the well of each contig
    """
    folder = str(Path(sample.output_dir).joinpath("temp_end_seqs"))
    os.makedirs(folder, exist_ok=True)

    # Index the contigs
    logging.info("Indexing assembled contigs... ")
    db = folder + "/contigs"
    blastdb_command = [executables["makeblastdb"]]
    blastdb_command += ["-in", sample.assembled_fosmids]
    blastdb_command += ["-dbtype", "nucl"]
    blastdb_command += ["-out", db]
    blastdb_command += ["1>/dev/null", "2>/dev/null"]
    subprocess_helper(blastdb_command)
    logging.info("done\n")

    logging.info("Aligning fosmid ends to assembly... ")
    blastn_command = [executables["blastn"]]
    blastn_command += ["-db", db]
    blastn_command += ["-outfmt", "\"6",
                       "sseqid", "slen", "qseqid", "pident", "length", "sstrand",
                       "sstart", "send", "bitscore", "qstart", "qend", "\""]
    blastn_command += ["-query", ends]
    blastn_command += ["-perc_identity", str(80)]
    blastn_command += ["-out", folder + f"/endsMapped.tbl"]
    blastn_command += ["1>/dev/null", "2>/dev/null"]
    subprocess_helper(blastn_command)
    logging.info("done\n")
    return


def get_fosmid_end_name(string):
    """
    This function is very unstable and currently needs to be manually changed for every FabFos run
    """
    i1, i2 = next(re.finditer(ENDS_NAME_REGEX, string)).span()
    name = string[i1:i2]
    # print(name)

    # print(string[i1:i2])
    # # PPSLIBM-08-H17_F.ab1
    # full_name = string.split('.')[0]
    # prefix = full_name.split('-')[0]
    # # print "string= ", string, "\tfull= ", full_name, "\tprefix= ", prefix
    # name, direction = full_name.split('_')
    # # name = prefix + '.' + name
    return name


def parse_end_alignments(sample: Sample, fosmid_assembly, ends_stats: FosmidEnds) -> dict:
    """
    Function parses BLAST alignments of fosmid ends mapped to the MEGAHIT assemblies

    :param sample: Miffed object with information of current sample
    :param fosmid_assembly: dictionary with headers as keys and sequences for values
    :param ends_stats: A FosmidEnds instance
    :return: dictionary of contig names (keys) and list of EndAlignment objects (values) or empty list
    """
    ends_mapping = dict()
    blast_output = str(Path(sample.output_dir).joinpath(f"temp_end_seqs/endsMapped.tbl"))
    try:
        blast_table = open(blast_output, 'r')
    except IOError:
        sys.exit("Unable to open " + blast_output)

    for contig in fosmid_assembly.keys():
        name = contig.split(" ")[0]
        ends_mapping[name] = list()
    line = blast_table.readline()
    n = 0
    while line:
        n += 1
        fields = line.strip().split("\t")
        alignment = EndAlignment(fields)
        contig = ">" + alignment.sseqid
        if alignment.length > 100 and alignment.bitscore > 500:
            alignment.parse_fosmid_end_name()
            if alignment.name not in ends_stats.failed:
                ends_stats.aligned.add(alignment.name)

            # Screen for the best alignment of a fosmid end based on bitscore
            try:
                hits = len(ends_mapping[contig])
            except KeyError:
                sys.exit("KeyError: " + contig + " not found in ends_mapping.keys()\n" + str(ends_mapping.keys()))
            if hits == 0:
                ends_mapping[contig].append(alignment)
            else:
                while hits > 0:
                    index = hits - 1
                    if ends_mapping[contig][index].qseqid == alignment.qseqid:
                        if ends_mapping[contig][index].bitscore < alignment.bitscore:
                            ends_mapping[contig].pop(index)
                    hits -= 1
                ends_mapping[contig].append(alignment)
        line = blast_table.readline()

    unaligned = ends_stats.all_clones.difference(ends_stats.aligned)
    for fosmid_end in unaligned:
        if fosmid_end not in ends_stats.failed:
            ends_stats.unaligned.add(fosmid_end)
    ends_stats.num_unaligned = len(ends_stats.unaligned)
    ends_stats.num_aligned = len(ends_stats.aligned)

    return ends_mapping


def reverse_complement(dna_sequence):
    rev_comp_seq = ""
    rev_seq = dna_sequence[::-1]
    for nuc in rev_seq:
        if nuc == "A" or nuc == "a":
            rev_comp_seq += "T"
        if nuc == "T" or nuc == "t":
            rev_comp_seq += "A"
        if nuc == "C" or nuc == "c":
            rev_comp_seq += "G"
        if nuc == "G" or nuc == "g":
            rev_comp_seq += "C"
        if nuc == "N" or nuc == "n":
            rev_comp_seq += "N"
    return rev_comp_seq


def find_mate_pairs(hit_list):
    pairs = list()
    direction_keeper = dict()
    for hit in hit_list:
        if hit.name in direction_keeper.keys() and hit.direction != direction_keeper[hit.name]:
            pairs.append(hit.name)
        else:
            direction_keeper[hit.name] = hit.direction
    return pairs


def get_exterior_positions(hits, clone):
    """
    Used to find the exterior positions of fosmid-end alignment on the contig.

    :param hits: List of EndAlignment objects
    :param clone: Name of the clone - presumably matches an EndAlignment.name
    :return: start and end positions
    """
    clone_alignments = 0
    positions_dict = dict()
    positions_dict["plus_strand"] = list()
    positions_dict["minus_strand"] = list()
    forward_aligned = False
    reverse_aligned = False
    contig = hits[0].sseqid
    for alignment in hits:
        # print alignment.get_info()
        if alignment.name == clone:
            clone_alignments += 1
            if alignment.direction == "F":
                forward_aligned = True
            if alignment.direction == "R":
                reverse_aligned = True
            if alignment.sstrand == "minus":
                positions_dict["minus_strand"].append(alignment.sstart)
                positions_dict["minus_strand"].append(alignment.send)
            if alignment.sstrand == "plus":
                positions_dict["plus_strand"].append(alignment.sstart)
                positions_dict["plus_strand"].append(alignment.send)

    # if not forward_aligned:
    #     logging.warning("Forward strand did not align in single-pair clone " + clone + "\n")
    # if not reverse_aligned:
    #     logging.warning("Reverse strand did not align in single-pair clone " + clone + "\n")
    # if len(positions_dict["minus_strand"]) == 0 or len(positions_dict["plus_strand"]) == 0:
    #     logging.warning("Forward and reverse mate-pairs of " + clone + " aligned to the same strand of " + contig +
    #                     "!\nThis sequence should not be trusted due to assembly or library preparation errors\n")

    positions = positions_dict["plus_strand"] + positions_dict["minus_strand"]
    start = min(positions)
    end = max(positions)
    return start, end

def assign_clones(ends_mapping, ends_stats: FosmidEnds, fosmid_assembly):
    """
    Assigns fosmid-end sequences to assembled contigs in fosmid_assembly

    :param ends_mapping:
    :param ends_stats: A FosmidEnds instance
    :param fosmid_assembly:
    :return:
    """
    clone_map = list()
    # multi_fosmid_map tracks all the fosmids each fosmid-end aligns to
    multi_fosmid_map = dict()
    unassigned_contigs = set()
    single_orphans = set()
    multi_pairs = set()
    single_pairs = set()
    multi_orphans = set()
    for contig in ends_mapping:
        hits = len(ends_mapping[contig])
        fosmid = FosmidClone(contig, fosmid_assembly[contig])
        if hits == 0:
            fosmid.clone = "None"
            fosmid.evidence = "No hit"
            fosmid.ends = "None"
            fosmid.strand = "None"
            clone_map.append(fosmid)
            unassigned_contigs.add(contig)
        if hits == 1:
            fosmid.clone = ends_mapping[contig][0].name
            fosmid.evidence = "single-orphan"
            fosmid.ends = ends_mapping[contig][0].direction
            fosmid.strand = ends_mapping[contig][0].sstrand
            if ends_mapping[contig][0].sstrand == "plus":
                fosmid.extract_sequence(ends_mapping[contig][0].sstart, len(fosmid.sequence))
                fosmid.sequence += 50 * 'N'
            if ends_mapping[contig][0].sstrand == "minus":
                fosmid.extract_sequence(ends_mapping[contig][0].sstart, 1)
                fosmid.sequence = (50 * 'N') + fosmid.sequence
            clone_map.append(fosmid)
            single_orphans.add(contig)
            if fosmid.clone not in multi_fosmid_map.keys():
                multi_fosmid_map[fosmid.clone] = list()
            multi_fosmid_map[fosmid.clone].append(fosmid)
        elif hits > 1:
            pairs = find_mate_pairs(ends_mapping[contig])
            clones_seen = list()
            while hits > 0:
                index = hits - 1
                hit = ends_mapping[contig][index]
                fosmid.clone = hit.name
                # Some fosmid-end alignments do not fit any of these categories
                if len(pairs) > 1 and fosmid.clone in pairs:
                    multi_pairs.add(contig)
                    fosmid.evidence = "multi-pair"
                    fosmid.ends = "Both"
                    fosmid.strand = hit.sstrand
                    start, end = get_exterior_positions(ends_mapping[contig], fosmid.clone)
                    fosmid.extract_sequence(start, end)
                if len(pairs) == 1 and fosmid.clone in pairs:
                    single_pairs.add(contig)
                    fosmid.evidence = "single-pair"
                    fosmid.ends = "Both"
                    fosmid.strand = hit.sstrand
                    # Extract sequence interior of the aligned mate-pairs
                    start, end = get_exterior_positions(ends_mapping[contig], fosmid.clone)
                    fosmid.extract_sequence(start, end)
                elif len(pairs) == 0:
                    multi_orphans.add(contig)
                    fosmid.evidence = "multi-orphan"
                    fosmid.ends = hit.direction
                    fosmid.strand = hit.sstrand
                    if hit.sstrand == "plus":
                        fosmid.extract_sequence(hit.sstart, len(fosmid.sequence))
                        fosmid.sequence += 50 * 'N'
                    elif hit.sstrand == "minus":
                        fosmid.extract_sequence(hit.sstart, 1)
                        fosmid.sequence = (50 * 'N') + fosmid.sequence
                    if fosmid.clone not in multi_fosmid_map.keys():
                        multi_fosmid_map[fosmid.clone] = list()
                    multi_fosmid_map[fosmid.clone].append(fosmid)
                elif len(pairs) != 0 and fosmid.clone not in pairs:
                    fosmid.ends = hit.direction
                    fosmid.evidence = "auxiliary-orphan-to-pair"
                    fosmid.strand = hit.sstrand
                    if hit.sstrand == "plus":
                        fosmid.extract_sequence(hit.sstart, len(fosmid.sequence))
                        fosmid.sequence += 50 * 'N'
                    elif hit.sstrand == "minus":
                        fosmid.extract_sequence(hit.sstart, 1)
                        fosmid.sequence = (50 * 'N') + fosmid.sequence
                # This avoids adding redundant single-pair and multi-pair FosmidClones to clone_map:
                if fosmid.evidence == "single-pair" or fosmid.evidence == "multi-pair":
                    if fosmid.clone not in clones_seen:
                        clones_seen.append(fosmid.clone)
                        clone_map.append(fosmid)
                else:
                    clones_seen.append(fosmid.clone)
                    clone_map.append(fosmid)

                fosmid = FosmidClone(contig, fosmid_assembly[contig])
                hits -= 1
    ends_stats.single_orphan = len(single_orphans)
    ends_stats.unassigned = len(unassigned_contigs)
    ends_stats.multi_pair = len(multi_pairs)
    ends_stats.single_pair = len(single_pairs)
    ends_stats.multi_orphan = len(multi_orphans)
    return clone_map, multi_fosmid_map


def prune_and_scaffold_fosmids(sample, clone_map, multi_fosmid_map):
    # fragments_tsv = sample.output_dir + os.sep + f"{sample.id}_potential_fosmid_fragments.tsv"

    # try:
    #     fragments = open(fragments_tsv, 'w')
    # except IOError:
    #     sys.exit("Unable to open " + fragments_tsv + " for writing!")

    for clone in multi_fosmid_map:
        if len(multi_fosmid_map[clone]) == 2:
            i = 0
            prefix = ""
            suffix = ""
            leading_node = ""
            trailing_node = ""
            direction = ""
            new_fosmid = None
            while i < len(clone_map):
                fosmid = clone_map[i]
                if fosmid.clone == clone:
                    # The new fosmid cannot be built from promiscuous ends. Immediately exit the loop in this case.
                    if fosmid.ends == direction:
                        new_fosmid = None
                        prefix = ""
                        suffix = ""
                        i = len(clone_map)
                    else:
                        direction = fosmid.ends
                        # Remove the old FosmidClone objects from clone_map
                        try:
                            mate = clone_map.pop(i)
                        except IndexError:
                            sys.exit("IndexError while parsing clone_map in prune_and_scaffold_fosmids:\n" +
                                     "pop index out of range: i = " + str(i) +
                                     ", len(clone_map) = " + str(len(clone_map)))
                        if mate.strand == "minus" and not suffix:
                            # Minus strand alignments have N's prepended
                            suffix = mate.sequence
                            trailing_node = mate.contig
                        if mate.strand == "minus" and suffix:
                            prefix = reverse_complement(mate.contig)
                            leading_node = mate.contig
                        if mate.strand == "plus" and not prefix:
                            # Plus strand alignments have N's appended
                            prefix = mate.sequence
                            leading_node = mate.contig
                        elif mate.strand == "plus" and prefix:
                            suffix = reverse_complement(mate.sequence)
                            trailing_node = mate.contig
                        # Scaffold the fosmids
                        if suffix and prefix:
                            new_fosmid = FosmidClone(leading_node + ',' + trailing_node, prefix + "-" + suffix)
                            new_fosmid.clone = clone
                            new_fosmid.evidence = "twin-orphan multi-contig"
                            new_fosmid.ends = "Both"
                            new_fosmid.strand = "NA"
                            new_fosmid.seq_length = len(prefix + suffix)
                            # Exit the loop
                            i = len(clone_map)
                i += 1
            if new_fosmid:
                clone_map.append(new_fosmid)
        else:
            # # Not confident in scaffolding these contigs
            # # Write the FosmidClone information to a potential_fosmid_fragments.tsv file
            # for fosmid in multi_fosmid_map[clone]:
            #     fragments.write(fosmid.get_info() + "\n")
            pass

    # fragments.close()

    return clone_map


def write_fosmid_assignments(sample, clone_map: list[FosmidClone], min_length: int):
    ws = Path(sample.output_dir)
    _fpaths = [ws.joinpath(f"contigs_{f}.fasta") for f in [
        "not_mapped", "single_mapped", "both_mapped", f"less_than_{min_length}",
    ]]
    _none, _single, _complete, _short = [open(p, "w") for p in _fpaths]

    def write_header(file, id: int, clone, l, mapping=None):
        m = f" {mapping}" if mapping is not None else ""
        c = f"_{clone}" if clone != "None" else ""
        file.write(f">{id:04}{c} length={l}{m}"+"\n")

    clone_map = sorted(clone_map, key=lambda c: c.seq_length, reverse=True)
    for i, fosmid_clone in enumerate(clone_map):
        id = i+1
        ends = {"F":"FW", "R":"RE"}.get(fosmid_clone.ends, fosmid_clone.ends)
        contig_seq = fosmid_clone.sequence

        if fosmid_clone.seq_length < min_length:
            write_header(_short, id, clone=str(fosmid_clone.clone), l=len(contig_seq), mapping=f"{fosmid_clone.evidence}|{ends}".replace(" ", "_").lower())
            _short.write(contig_seq)
            _short.write("\n")
        else:
            # write_header(_all, id, clone=str(fosmid_clone.clone), l=len(contig_seq), mapping=f"{fosmid_clone.evidence}|{ends}".replace(" ", "_").lower())
            fasta = {
                "None": _none,
                "Both": _complete,
                "FW": _single,
                "RE": _single,
            }[ends]
            write_header(fasta, id, clone=str(fosmid_clone.clone), l=len(contig_seq), mapping=ends if ends in {"FW", "RE"} else None)

            # for fa in [_all, fasta]:
            for fa in [fasta]:
                fa.write(contig_seq)
                fa.write("\n")

    # for file_handle in _all, _none, _single, _complete:
    for file_handle in _none, _single, _complete:
        file_handle.close()
    return _fpaths[0]

def read_fastq(file_path):
    fq = SeqIO.parse(file_path, "fastq")
    for entry in fq:
        yield entry.id, entry.seq

def read_fasta(fasta):
    """
    Function to read a fasta file and store it in a dictionary. Headers are keys, sequences are values

    :param fasta: A fasta file name (and path, if necessary)
    :return: dictionary with headers as keys and sequences for values
    """
    # TODO: Replace this with the much faster pyfastx API
    fasta_dict = dict()
    header = ""
    sequence = ""
    fasta_seqs = open(fasta, 'r')
    for line in fasta_seqs:
        line = line.replace("\n", "").strip()
        if len(line) == 0: continue
        if line[0] == '>':
            if header != "":
                fasta_dict[header] = sequence
            # Takes the first word in the header (to match alignment outputs)
            header = line.split(" ")[0]
            sequence = ""
        else:
            sequence += line
    fasta_dict[header] = sequence

    fasta_seqs.close()

    return fasta_dict


def get_assembler_version(assembler):
    """
    Wrapper function for retrieving the version number of either SPAdes or MEGAHIT

    :return: version string
    """
    if assembler == "spades":
        ver_proc = subprocess.Popen(["spades.py", "-v"], stdout=subprocess.PIPE)
    elif assembler == "megahit":
        ver_proc = subprocess.Popen(["megahit", "-v"], stdout=subprocess.PIPE)
    else:
        logging.error("Version cannot currently be identified for assembler '" + assembler + "'.\n")
        sys.exit(3)
    out, err = ver_proc.communicate()
    asm_ver = str(out.strip(), encoding="utf-8")
    return asm_ver


def get_fasta_size(fasta_dict: dict) -> int:
    return sum([len(seq) for seq in fasta_dict.values()])


def get_assembly_nx(fasta_dict: dict, increment=0.1) -> dict:
    """
    Calculates Nx stats of an assembly stored in a dictionary according to some incrementing float

    :param fasta_dict: A dictionary containing sequence names as keys and sequences as values
    :param increment: A float to increment the proportion of the assembly to be calculated
    :return: A dictionary mapping proportions (0-1) to the smallest sequence length required to sum to that proportion
    """
    nx_stats = {}
    asm_len = get_fasta_size(fasta_dict)
    proportion = 0.0
    running_sum = 0
    for name in reversed(sorted(fasta_dict, key=lambda x: len(fasta_dict[x]))):
        seq = fasta_dict[name]
        running_sum += len(seq)
        while running_sum >= (asm_len*proportion) and proportion <= 1:
            nx_stats[round(proportion, 1)] = len(seq)
            proportion += increment
    return nx_stats


def get_assembly_stats(fosmid_assembly: dict, assembler: str, ):
    assembly_stats = dict()
    assembly_stats["Assembler"] = get_assembler_version(assembler)

    nx_stats = get_assembly_nx(fosmid_assembly)

    # Print some stats to the log
    logging.info("N50 = " + str(nx_stats[0.5]) + "\n")
    logging.info("Longest contig = " + str(nx_stats[0]) + "bp\n")
    # Save Nx stats for the output
    assembly_stats["N50"] = str(nx_stats[0.5])
    assembly_stats["N90"] = str(nx_stats[0.9])

    size = 0
    num_27k = 0
    num_50k = 0
    for contig in fosmid_assembly:
        sequence = fosmid_assembly[contig]
        seq_len = len(sequence)
        size += seq_len
        if seq_len >= 27000:
            num_27k += 1
            if seq_len >= 50000:
                num_50k += 1

    assembly_stats["27kbp"] = str(num_27k)
    assembly_stats["50kbp"] = str(num_50k)
    assembly_stats["Contigs"] = str(len(fosmid_assembly))
    assembly_stats["basepairs"] = str(size)

    return assembly_stats


# def update_metadata(metadata_file, library_metadata: Miffed,
#                     read_stats: ReadStats, assembly_stats, ends_stats: FosmidEnds):
#     try:
#         metadata = open(metadata_file, 'a')
#     except IOError:
#         logging.error("Unable to open " + metadata_file + " to append library metadata!\n")
#         sys.exit(3)

#     metadata.write(library_metadata.id + "\t")
#     metadata.write(library_metadata.project + "\t")
#     metadata.write(library_metadata.selector + "\t")
#     metadata.write(library_metadata.vector + "\t")
#     metadata.write(library_metadata.screen + "\t")
#     metadata.write(library_metadata.selection + "\t")
#     metadata.write(str(library_metadata.num_fosmids_estimate) + "\t")
#     metadata.write(library_metadata.seq_submission_date + "\t")
#     metadata.write(library_metadata.glycerol_plate + "\t")
#     metadata.write(library_metadata.seq_center + "\t")
#     metadata.write(library_metadata.seq_type + "\t")
#     metadata.write(str(library_metadata.read_length) + "\t")
#     metadata.write(library_metadata.instrument + "\t")
#     metadata.write(strftime("%Y-%m-%d") + "\t")
#     metadata.write(str(read_stats.num_raw_reads) + "\t")
#     metadata.write(str(read_stats.percent_on_target()) + "\t")
#     metadata.write(str(read_stats.num_reads_trimmed) + "\t")
#     metadata.write(assembly_stats["Assembler"] + "\t")
#     metadata.write(assembly_stats["Contigs"] + "\t")
#     metadata.write(assembly_stats["N50"] + "\t")
#     metadata.write(assembly_stats["N90"] + "\t")
#     metadata.write(assembly_stats["27kbp"] + "\t")
#     metadata.write(assembly_stats["50kbp"] + "\t")
#     metadata.write(str(ends_stats.single_pair) + "\t")
#     metadata.write(str(ends_stats.multi_pair) + "\t")
#     metadata.write(str(ends_stats.single_orphan) + "\t")
#     metadata.write(str(ends_stats.multi_orphan) + "\t")
#     metadata.write(str(ends_stats.unassigned) + "\t")
#     metadata.write(str(ends_stats.num_total) + "\t")
#     metadata.write(str(ends_stats.num_aligned) + "\t")
#     metadata.write(str(ends_stats.num_unaligned) + "\t")
#     metadata.write(str(ends_stats.num_failed) + "\n")

#     metadata.close()
#     return


def write_fosmid_end_failures(sample: Sample, ends_stats: FosmidEnds) -> None:
    """
    Writes the names of fosmid ends that could not be aligned and failed

    :param sample: Miffed object with information of current sample
    :param ends_stats: An instance of the FosmidEnds class
    :return:
    """
    sample_prefix = sample.output_dir + os.sep + sample.id
    fosmid_end_failures = sample_prefix + "_end_mapping_failures.tsv"

    try:
        failure_file = open(fosmid_end_failures, 'w')
    except IOError:
        logging.error("Cannot open file: " + fosmid_end_failures)
        sys.exit(3)
    failure_file.write("#Fosmid-end\tcategory\n")
    for end in ends_stats.failed:
        failure_file.write(end + "\tSequencing failure\n")
    for end in ends_stats.unaligned:
        failure_file.write(end + "\tUnaligned\n")

    failure_file.close()

    # For comparing the fosmid-end names in each set:
    # i = 0
    # while i < 5:
    #     print "all", ends_stats["all_clones"].pop()
    #     print "failed", ends_stats["failed"].pop()
    #     print "unaligned", ends_stats["unaligned"].pop()
    #     print "aligned", ends_stats["aligned"].pop()
    #     i += 1
    #
    # with open("tmp", 'w') as tmp:
    #     for x in ends_stats["failed"]:
    #         tmp.write("failed " + x + "\n")
    #     for x in ends_stats["unaligned"]:
    #         tmp.write("unaligned " + x + "\n")
    #     for x in ends_stats["aligned"]:
    #         tmp.write("aligned " + x + "\n")

    return


# def write_unique_fosmid_ends_to_bulk(args):
#     logging.info("Writing new fosmid end sequences to the legacy fosmid end file... ")
#     new_fosmid_ends = read_fasta(args.ends)
#     bulk_ends_file = args.output + "/FabFos_legacy_ends.fasta"
#     bulk_ends = read_fasta(bulk_ends_file)

#     try:
#         legacy_ends = open(bulk_ends_file, 'a')
#     except IOError:
#         raise IOError("Cannot open " + bulk_ends_file + " for appending!")

#     for new_end in new_fosmid_ends:
#         if new_end not in bulk_ends.keys():
#             if len(new_fosmid_ends[new_end]) > 100:
#                 legacy_ends.write(new_end + "\n")
#                 legacy_ends.write(new_fosmid_ends[new_end] + "\n")

#     legacy_ends.close()
#     logging.info("done.\n")

#     return


def align_nanopore_to_background(args, executables, sample):
    logging.info("Aligning " + sample.nanopore + " to background sequences... ")

    folder = str(Path(sample.output_dir).joinpath("temp_nanopore"))
    # Make the BLAST database
    try:
        background_blastdb_stdout = open(folder + os.sep + "background_blastdb.stdout", 'w')
    except IOError:
        logging.error("Unable to open " + folder + os.sep + "background_blastdb.stdout for writing")
        sys.exit(3)

    background_db = args.background + "_BLAST"
    blast_db_command = [executables["makeblastdb"]]
    blast_db_command += ["-dbtype", "nucl"]
    blast_db_command += ["-in", args.background]
    blast_db_command += ["-out", background_db]

    stdout, retcode = subprocess_helper(cmd_list=blast_db_command)

    background_blastdb_stdout.write(stdout)
    background_blastdb_stdout.close()

    # Align the corrected reads to the BLAST database
    nanopore_background_alignments = folder + os.sep + "nanopore_background_BLAST.tbl"
    blastn_command = [executables["blastn"]]
    blastn_command += ["-query", sample.nanopore]
    blastn_command += ["-db", background_db]
    blastn_command += ["-outfmt", "\"6",
                       "sseqid", "slen", "qseqid", "pident", "length", "sstrand",
                       "sstart", "send", "bitscore", "qstart", "qend", "\""]
    blastn_command += ["-out", nanopore_background_alignments]
    subprocess_helper(blastn_command)

    # Remove the BLAST database
    blastdb_files = glob.glob('.'.join(args.background.split('.')[0:-1]) + "_BLAST.*")
    for db_file in blastdb_files:
        os.remove(db_file)

    logging.info("done.\n")
    return nanopore_background_alignments


def filter_background_nanopore(sample: Sample, nanopore_background_alignments) -> dict:
    """
    Function for removing background (i.e. vector backbone, genome) reads specifically for Nanopore long reads

    :param nanopore_background_alignments:
    :param sample: Miffed object with information of current sample
    :return: nanopore_reads is a dictionary of trimmed nanopore sequence reads
    """

    nanopore_reads = read_fasta(sample.nanopore)
    logging.info("Nanopore reads prior to trimming background sequence = " +
                 str(len(nanopore_reads.keys())) + "\n")
    reads_with_background = dict()
    truncated = 0
    if sample.error_correction is False:
        identity_threshold = 85
    else:
        identity_threshold = 95

    try:
        alignments = open(nanopore_background_alignments, 'r')
    except IOError:
        raise IOError("Unable to open " + nanopore_background_alignments + " for reading!")

    for line in alignments:
        fields = line.strip().split('\t')
        hit = Alignment(fields)
        # Only working with high confidence alignments
        if hit.pident < identity_threshold:
            continue
        if '>' + hit.qseqid not in reads_with_background.keys():
            reads_with_background['>' + hit.qseqid] = list()
        reads_with_background['>' + hit.qseqid].append(hit)
        # Skip the reads where the hit.length + 1000 is longer than the read length
        if hit.length + 1000 > len(nanopore_reads['>' + hit.qseqid]):
            nanopore_reads['>' + hit.qseqid] = ""

    # Find the most interior positions that aligned to background sequence with high identity
    for background_read in reads_with_background:
        read_seq = nanopore_reads[background_read]
        if len(read_seq) > 0:
            truncated += 1
            five_prime, three_prime = get_internal_positions(reads_with_background[background_read])
            # Trim
            if three_prime > 1000:
                nanopore_reads[background_read] = read_seq[:three_prime]
            if three_prime <= 1000:
                nanopore_reads[background_read] = read_seq[five_prime:]

    logging.info("Nanopore reads that were trimmed for background sequence = " + str(truncated) + "\n")
    sample.read_stats.num_on_target = truncated
    sample.read_stats.num_reads_trimmed = truncated

    return nanopore_reads


def get_internal_positions(alignments):
    """
    Function for returning the most internal positions from a list of alignments

    :param alignments:
    :return:
    """
    five_prime = 0
    three_prime = 0
    for hit in alignments:
        if hit.qstart < hit.qend:
            if hit.qend > five_prime:
                five_prime = hit.qend
                three_prime = hit.qstart
        if hit.qend < hit.qstart:
            raise AssertionError(" ".join(["Unexpected 3'-5' alignment by BLAST:",
                                           hit.qseqid, hit.sseqid, hit.qstart, hit.qend]))
    return five_prime, three_prime


def write_trimmed_reads(miffed_entry, nanopore_reads):
    majority_background = 0
    trimmed_nanopore_fasta = str(Path(miffed_entry.output_dir).joinpath("nanopore_trimmed.fasta"))
    try:
        new_fasta = open(trimmed_nanopore_fasta, 'w')
    except IOError:
        logging.error("Unable to open " + trimmed_nanopore_fasta + " for writing.\n")
        sys.exit()

    for trimmed_read in nanopore_reads:
        if len(nanopore_reads[trimmed_read]) >= 1000:
            new_fasta.write(trimmed_read + "\n")
            new_fasta.write(nanopore_reads[trimmed_read] + "\n")
        else:
            majority_background += 1

    new_fasta.close()

    logging.info("Nanopore reads that are shorter than 1000bp after trimming background = " +
             str(majority_background) + "\n")

    return


def determine_min_count(num_reads: int, num_fosmids: int, k_max: int):
    """
    Function to determine the best value for the minimum coverage of a k-mer to be included in assembly

    :param num_reads:
    :param num_fosmids:
    :param k_max:
    :return:
    """
    approx_coverage = (num_reads * k_max) / (int(num_fosmids) * 40000)
    sys.stdout.write("Estimated fosmid coverage = " + str(approx_coverage) + "\n")
    sys.stdout.flush()
    min_count = 10
    dist_tail = approx_coverage / 100
    if dist_tail > min_count:
        min_count = int(dist_tail)

    return min_count, approx_coverage


def assemble_nanopore_reads(sample: Sample, canu_exe, skip_correction, threads=2):
    """
    Wrapper function for launching canu assembler with the corrected nanopore reads

    :param sample:
    :param canu_exe:
    :param skip_correction:
    :param threads:
    :return:
    """
    logging.info("Assembling error-corrected nanopore reads with canu... ", "out")

    canu_output = sample.output_dir + os.sep + "canu" + os.sep
    os.makedirs(canu_output, exist_ok=True)

    corrected_nanopore_fasta = sample.output_dir + os.sep + sample.id + "_nanopore_trimmed.fasta"
    canu_stdout = open(canu_output + "canu.stdout", 'w')

    genome_size = int(sample.num_fosmids_estimate) * 40

    canu_command = [canu_exe]
    canu_command += ["genomeSize=" + str(genome_size) + "k", "minThreads=" + str(threads)]
    if sample.error_correction is True:
        # TODO: Compare assemblies where proovread-corrected nanopore reads are just assembled or ran through all stages
        canu_command.append("-assemble")
    canu_command += ["-p", sample.id]
    canu_command += ["-d", canu_output]
    if skip_correction:
        canu_command += ["-nanopore-raw", corrected_nanopore_fasta]
    else:
        canu_command += ["-nanopore-corrected", corrected_nanopore_fasta]

    stdout, retcode = subprocess_helper(cmd_list=canu_command)
    canu_stdout.write(stdout)
    canu_stdout.close()

    sample.assembled_fosmids = str(Path(canu_output).joinpath(f"{sample.id }.contigs.fasta"))
    logging.info("done.\n")
    return

def fabfos_main(sys_args):
    args = get_options(sys_args)

    executables = find_executables(assembler=args.assembler, nanopore=args.nanopore_reads)
    versions_dict = find_dependency_versions(executables)
    validate_dependency_versions(versions_dict)
    summarize_dependency_versions(dep_versions=versions_dict)
    if args.version:
        sys.exit(0)

    fos_father = FabFos(args.output)
    review_arguments(args, fos_father)
    
    out_path = Path(args.output)
    if out_path.exists():
        if args.overwrite:
            # logging.info(f"cleaning output folder [{out_path}]\n")
            # time.sleep(3)
            os.system(f"rm -r {out_path}")
        else:
            logging.error(f"output folder [{out_path}] already exists\n")
            sys.exit(1)

    # create temp workspace
    fos_father.furnish()
    # Setup the global logger and main log file
    prep_logging(args.output + os.sep + "fabfos.log", args.verbose)

    ends_stats = FosmidEnds()
    if args.ends:
        ends_stats.fasta_path = args.ends
        ends_stats.load_ends()

    sample_id = out_path.name
    sample = Sample(sample_id)
    ghetto_sync_args(args, fos_father, sample)

    # read qc
    sample.gather_reads(args.reads, args.reverse, args.parity, executables, sample.output_dir, args.type)
    Fastqc(out_path, [Path(r) for r in [sample.forward_reads, sample.reverse_reads] if r is not None])
    qc_pe, qc_se = sample.qc_reads(args.background, args.parity, fos_father.adaptor_trim, executables, args.threads)
    qced_reads = qc_pe if len(qc_pe) > 0 else [qc_se[0]]
    
    # pool size estimate
    if args.pool_size is None and args.vector:
        logging.info(f"Estimating fosmid pool size... ")
        size_estimate = EstimateFosmidPoolSize([Path(r) for r in qced_reads], Path(args.vector), Path(args.output).joinpath("temp_pool_size_estimate"), args.threads)
        if size_estimate is None:
            logging.error(f"failed to estimate pool size, defaulting to given estimate\n")
            if args.pool_size is None:
                logging.error(f"no given estimate\n")
                sys.exit(1)
            else:
                size_estimate = args.pool_size
        else:
            logging.info(f"done.\n")
        logging.info(f"Pool size estimate = {size_estimate}\n")
    else:
        size_estimate = args.pool_size
    sample.num_fosmids_estimate = size_estimate # used in determine_min_count() of assemble_fosmids()
    
    if sample.nanopore:
        sample.prep_nanopore(args, executables, args.nanopore_reads)
        # Assemble nanopore reads
        assemble_nanopore_reads(sample, executables["canu"], args.skip_correction, args.threads)
    else:
        sample.assemble_fosmids(executables, args.threads)

    if not Path(sample.assembled_fosmids).exists():
        logging.error(f"assembly failed")
        sys.exit(1)

    # trim vector
    if args.vector:
        sample.assembled_fosmids = TrimBackbone(out_path, Path(args.vector), Path(sample.assembled_fosmids))

    fosmid_assembly = read_fasta(sample.assembled_fosmids)
    assembly_stats = get_assembly_stats(fosmid_assembly, sample.assembler)
    assembly_stats["est_n_fosmids"] = size_estimate

    # For mapping fosmid ends:
    if args.ends:
        logging.info(f"Mapping end sequences:\n")
        map_ends(executables, args.ends, sample)
        ends_mapping = parse_end_alignments(sample, fosmid_assembly, ends_stats)
        clone_map, multi_fosmid_map = assign_clones(ends_mapping, ends_stats, fosmid_assembly)
        clone_map = prune_and_scaffold_fosmids(sample, clone_map, multi_fosmid_map)
        final_contig_file = write_fosmid_assignments(sample, clone_map, int(args.min_length))
        write_fosmid_end_failures(sample, ends_stats)
    else:
        # some function above does this if ends are provided and also adds end mappings to the headers
        sample.assembled_fosmids = str(FilterMinLength(Path(sample.assembled_fosmids), int(args.min_length)))
        final_contig_file = Path(sample.output_dir).joinpath(f"{sample.id}_contigs.fasta")
        shutil.copy(sample.assembled_fosmids, final_contig_file)

    is_pe = args.parity != "se"
    reads = qc_pe if is_pe else qc_se
    AssemblyStats(out_path, [Path(r) for r in reads], final_contig_file, int(args.threads), paired_end=is_pe)

    # clean up, organize outputs
    logging.info(f"Final cleanup... ")
    ws = Path(sample.output_dir)
    meta = {}

    # - read stats
    tmp_cov = ws.joinpath("estimated_coverage.txt")
    if tmp_cov.exists():
        with open(tmp_cov) as f:
            cov = f.readline()[:-1]
            meta["est_fosmid_coverage"] = cov
    meta["n_reads"] = sample.read_stats.num_raw_reads
    meta["n_reads_after_qc"] = sample.read_stats.num_filtered_reads
    os.unlink(tmp_cov)

    # # - asm stats
    # def _parse_k(k):
    #     if k.endswith("kbp"): k = f"atleast_{k}"
    #     return k.lower()
    # meta |= {f"asm_{_parse_k(k)}":v for k, v in assembly_stats.items()}
    
    # # - write meta to table
    # with open(ws.joinpath(f"{sample.id}_metadata.tsv"), "w") as f:
    #     header, values = [], []
    #     for k, v in meta.items():
    #         header.append(str(k))
    #         values.append(str(v))
    #     for l in [header, values]:
    #         f.write("\t".join(l)+"\n")
    os.system(f"""\
        mv {ws}/fastqc {ws}/{sample_id}_fastqc
        mv {ws}/quast {ws}/{sample_id}_quast
        mv {ws}/assembly_coverage.tsv {ws}/{sample_id}_read_coverage.tsv
    """)

    # - write end mapping to table
    end_map = ws.joinpath("end_seqs/endsMapped.tbl")
    if end_map.exists():
        with open(ws.joinpath(f"{sample_id}_end_mapping.tsv"), "w") as out:
            HEADER = [
                "sseqid", "slen", "qseqid", "pident", "length", "sstrand",
                "sstart", "send", "bitscore", "qstart", "qend"
            ]
            out.write("\t".join(HEADER)+"\n")
            with open(end_map) as f:
                for l in f: out.write(l)

    logging.info(f"done.\n")
    logging.info(f"Results are in: {sample.output_dir}\n")
