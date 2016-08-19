#!/usr/bin/env python

try:
    import argparse
    import sys
    import glob
    import os
    import traceback
    import re
    import subprocess
    import shutil
    from time import strftime
except ImportWarning:
    sys.stderr.write("Could not load some user defined module functions")
    sys.stderr.write(traceback.print_exc(10))
    sys.exit(3)

"""
FabFos: a pipeline for automatically performing quality controls, assembly, and data storage management
for fosmid sequence information. Circa 2016 - Hallam Lab, UBC
"""

__version__ = "0.1"
__author__ = "Connor Morgan-Lang"
__license__ = "GPL"
__maintainer__ = "Connor Morgan-Lang"
__email__ = "c.morganlang@gmail.com"
__status__ = "Unstable"


def get_dict_keys(dictionary):
    """
    A function for returning the sorted list of dictionary keys, regardless of the Python version.
    :param dictionary: dict() object
    :return: list() of sorted dictionary keys
    """
    dictionary.keys()
    return


def stdprint(obj, channel, cap=""):
    str_obj = str(obj)
    if type(str_obj) is not str:
        sys.exit("TypeError: stdprint only accepts string objects!")
    elif str(channel) == "err":
        sys.stderr.write(str_obj + cap)
        sys.stderr.flush()
    elif str(channel) == "out":
        sys.stdout.write(str_obj + cap)
        sys.stdout.flush()
    else:
        sys.stderr.write("ERROR: Unrecognized input to stdprint:\n")
        sys.stderr.write(str_obj + "\n" + channel + "\n")
        sys.exit()


def get_options():
    parser = argparse.ArgumentParser(description="Pipeline for processing and organizing fosmid seqeunce information."
                                                 "NOTE: maximum of 2 sequence files--user must concatenate files if >2")
    parser.add_argument("-m", "--miffed", type=str, required=True,
                        help="The minimum information for fosmid environmental data (e.g., sample ID, "
                             "sequencing platform, environment, project) in a comma-separated file. [.csv]")
    parser.add_argument("-b", "--backbone", type=str, required=True,
                        help="Path to the fosmid backbone database [.fasta]")
    parser.add_argument("-r", "--reads", type=str, required=True,
                        help="Path to the sequenced reads directory. This parameter is for the single-end reads "
                             "(Sanger), the forward strand file, or the interleaved paired-end file. [.fastq]")
    parser.add_argument("-f", "--fabfos_path", type=str, required=False,
                        default="/mnt/nfs/sharknado/LimsData/FabFos/",
                        help="Path to FabFos database on sharknado [DEFAULT = /mnt/nfs/sharknado/LimsData/FabFos/]")
    # TODO: Determine how to handle interleaved paired-end files
    parser.add_argument("-T", "--threads", type=str, required=False, default=str(8),
                        help="The number of threads that can be used [DEFAULT = 8]")
    parser.add_argument("-2", "--reverse", type=str, required=False,
                        help="Path to the directory containing reverse-end reads (if applicable) [.fastq]")
    parser.add_argument("-i", "--interleaved", required=False, default=False, action="store_true",
                        help="Flag indicating the reads are interleaved "
                             "(i.e. forward and reverse pairs are in the same file)")
    parser.add_argument("-e", "--ends", required=False, default=None,
                        help="FASTA file containing fosmid ends - these will be used for alignment to fosmid contigs."
                             " [.fasta]")
    args = parser.parse_args()
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


def review_arguments(args):
    """
    Function that ensures the proper information has been provided and adds system information to args
    :param args: parsed command-line arguments from get_options()
    :return: args with more information
    """
    # Review the provided arguments:
    if not os.path.isdir(args.reads):
        raise IOError("ERROR: " + args.reads + " is not a valid directory!")
    if args.reverse and not os.path.isdir(args.reverse):
        raise IOError("ERROR: " + args.reverse + " is not a valid directory!")
        # raise FileNotFoundError("ERROR: " + args.reads + " does not exist!")

    if not os.path.isfile(args.backbone):
        sys.exit("ERROR: " + args.backbone + " does not exist!")
    # Add new information:
    args.os = os_type()
    if sys.version_info > (2, 9):
        args.py_version = 3
    else:
        args.py_version = 2

    args.adaptor_trim = "/home/cmorganlang/bin/Trimmomatic-0.35/adapters/"

    args.metadata_header = "#Sample Name (LLLLL-PP-WWW)\tProject\tHuman selector\tVector Name\t" \
                           "Screen [in silico | functional]\tSelection criteria\tNumber of fosmids\t" \
                           "Sequencing submission date (YYYY-MM-DD)\tGlycerol plate name\tSequencing center\t" \
                           "Sequencing type\tRead length\tInstrument\tDate of FabFos analysis\tNumber of reads\t" \
                           "% off-target reads\tNumber of trimmed reads\tNumber of Contigs\tN50\tN90\t" \
                           "Contigs > 27kbp\tContigs > 50kbp\t" \
                           "# single-pair\t# multi-pair\t# single-orphan\t# multi-orphan\t# Unidentifiable\t" \
                           "# Fosmid ends(X2 = #reads)\t# aligned fosmid ends\t# unaligned fosmid ends\t# Failed ends\n"

    if not os.path.isfile(args.fabfos_path + os.sep + "FabFos_master_metadata.tsv"):
        sys.exit("ERROR: This fabfos_path directory does not contain a FabFos_master_metadata.tsv file."
                 "Are you sure its a bona fide FabFos repository?")

    args.master_metadata = args.fabfos_path + os.sep + "FabFos_master_metadata.tsv"

    return args


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
            if is_exe(exe_file):
                return exe_file
    return None


def find_executables(args):
    """
    Find the executables that are to be used in the pipeline:
    BWA, MEGAHIT, Trimmomatic, samtools
    :param args: parsed command-line arguments from get_options()
    :return: args with executables list included
    """
    required_execs = ["bwa", "samtools", "megahit", "bam2fastq", "trimmomatic", "fq2fa", "blastn", "makeblastdb",
                      "splitFASTA", "getNx"]
    args.executables = dict()
    for executable in required_execs:
        args.executables[executable] = which(executable)
        if args.executables[executable] is None:
            raise EnvironmentError("ERROR: unable to find executable for " + executable)
    return args


def check_index(path, bwa_path):
    extensions = ['.bwt', '.pac', '.ann', '.amb', '.sa']
    index_parts = [path+e for e in extensions]
    for part in index_parts:
        if not os.path.isfile(part):
            stdprint("\nUnable to find indexed backbone database. Indexing now... ", "err")
            index_command = [bwa_path, "index", path]
            index_dir = os.path.dirname(os.path.abspath(path))
            index_command += ["1>", index_dir + os.sep + "bwa_index.stdout"]
            index_command += ["2>", index_dir + os.sep + "bwa_index.stderr"]
            p_index = subprocess.Popen(' '.join(index_command), shell=True, preexec_fn=os.setsid)
            p_index.wait()
            if p_index.returncode != 0:
                sys.exit("ERROR: bwa index did not complete successfully")
            stdprint("done.\nResuming alignment... ", "err")
            break
    return


def filter_backbone(sample_id, args, raw_reads):
    """
    Function to generate fastq files that do not contain any sequences in `backbone`.
    Depends on bwa, samtools, and bam2fastq
    :param sample_id: The unique identifier for this fosmid library
    :param args: parsed command-line arguments from get_options()
    :param raw_reads: Dictionary containing forward and reverse fastq file names
    :return: list of fastq files containing the filtered reads
    """
    stdprint("Filtering off-target reads... ", "out")
    check_index(args.backbone, args.executables["bwa"])
    align_command = [args.executables["bwa"], "mem", "-t", args.threads, args.backbone, raw_reads["forward"]]
    if args.reverse:
        align_command.append(raw_reads["reverse"])
    align_command.append("1>")
    sam_file = args.output_dir[sample_id] + sample_id + ".sam"
    align_command.append(sam_file)
    align_command += ["2>", "/dev/null"]
    p_mem = subprocess.Popen(' '.join(align_command), shell=True, preexec_fn=os.setsid)
    p_mem.wait()

    # Use samtools to convert sam to bam
    bam_convert = [args.executables["samtools"], "view", "-bS", "-@", args.threads, sam_file, "|"]
    # Piping to samtools sort
    bam_file = args.output_dir[sample_id] + os.sep + sample_id + "_sorted"
    bam_convert.append(args.executables["samtools"])
    bam_convert += ["sort", "-@", args.threads, "-", bam_file]
    p_samtools_stdout = open(args.output_dir[sample_id] + os.sep + "samtools.stdout", 'w')
    p_samtools_stderr = open(args.output_dir[sample_id] + os.sep + "samtools.stderr", 'w')
    p_samtools = subprocess.Popen(' '.join(bam_convert), shell=True, preexec_fn=os.setsid,
                                  stdout=p_samtools_stdout, stderr=p_samtools_stderr)
    p_samtools.wait()
    p_samtools_stdout.close()
    p_samtools_stderr.close()

    # extract the unaligned reads using bam2fastq
    fastq_extract = [args.executables["bam2fastq"], "-o",
                     args.output_dir[sample_id] + os.sep + sample_id + "_BackboneFiltered_R#.fastq", "--no-aligned",
                     bam_file + ".bam"]
    fastq_extract += ["1>", args.output_dir[sample_id] + os.sep + "bam2fastq.stdout"]
    fastq_extract += ["2>", args.output_dir[sample_id] + os.sep + "bam2fastq.stderr"]
    p_bam2fastq = subprocess.Popen(' '.join(fastq_extract), shell=True, preexec_fn=os.setsid)
    p_bam2fastq.wait()
    filtered_reads = glob.glob(args.output_dir[sample_id] + os.sep + sample_id + "_BackboneFiltered*fastq")
    stdprint("done.", "out", "\n")
    return filtered_reads


def quality_trimming(sample_id, args, filtered_reads):
    """
    Wrapper for trimmomatic
    :param sample_id: The unique identifier for the sample - used for output directory creation and file naming
    :param args: command-line arguments list
    :param filtered_reads: list of backbone-filtered fastq files
    :return: list of quality-trimmed fastq files
    """
    stdprint("Trimming reads... ", "out")
    trimmomatic_command = ["java", "-jar", args.executables["trimmomatic"], "PE", "-threads", args.threads]
    trimmomatic_command += filtered_reads
    trim_prefix = args.output_dir[sample_id] + os.sep + sample_id + "_trim_"
    trimmomatic_outputs = [trim_prefix + "pe.1.fq", trim_prefix + "se.1.fq",
                           trim_prefix + "pe.2.fq", trim_prefix + "se.2.fq"]
    trimmomatic_command += trimmomatic_outputs
    # (/home/cmorganlang/bin/Trimmomatic-0.35/adapters/TruSeq3-PE.fa)
    adapters = "ILLUMINACLIP:" + args.adaptor_trim + "TruSeq3-PE.fa:2:3:10"
    trimmomatic_command.append(adapters)
    trimmomatic_command += ["LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"]
    try:
        trim_stderr = open(trim_prefix+"stderr.txt", 'w')
    except:
        raise IOError("ERROR: cannot open file: " + trim_prefix + "stderr.txt")
    try:
        trim_stdout = open(trim_prefix + "stdout.txt", 'w')
    except:
        raise IOError("ERROR: cannot open file: " + trim_prefix + "stdout.txt")
    p_trim = subprocess.Popen(' '.join(trimmomatic_command), shell=True, preexec_fn=os.setsid,
                              stderr=trim_stderr, stdout=trim_stdout)
    p_trim.wait()
    trim_stderr.close()
    trim_stdout.close()
    stdprint("done.", "out", "\n")
    return trimmomatic_outputs


def assemble_fosmids(sample_id, args, reads, k_min, k_max,  trimmed_reads):
    """
    Wrapper function for MEGAHIT assembler - multi-sized de Bruijn graph based assembler for metagenomes
    :param sample_id: The unique identifier for the sample - used for output directory creation
    :param args: command-line arguments list
    :param reads: List containing paired-end and single-end input reads for MEGAHIT
    :param k_min: 71; painstakingly determined to be optimal for fosmids sequenced with Illumina
    :param k_max: 241; painstakingly determined to be optimal for fosmids sequenced with Illumina
    :param trimmed_reads: list of files output by trimmomatic (returned by quality_trimming)
    :return:
    """
    stdprint("Assembling sequences using MEGAHIT.", "out", "\n")
    stdprint("Parameters:\n--k-min = " + str(k_min) + "\t--k-max = " + str(k_max) + "\t--k-step = 10\t--min-count = 10",
             "out",
             "\n")
    se = ""
    pe = ""
    for fasta in reads:
        if re.search(r"_paired.fa", fasta):
            pe = fasta
        else:
            se = fasta
    for fastq in trimmed_reads:
        if re.search(r'pe.1.fq$', fastq):
            forward = fastq
        if re.search(r'pe.2.fq$', fastq):
            reverse = fastq
    if not os.path.isfile(se):
        sys.exit("ERROR: FASTA file containing trimmed orphaned reads could not be located!")
    if not os.path.isfile(pe):
        sys.exit("ERROR: FASTA file containing trimmed paired-end reads could not be located!")

    # spades_command = ["/home/cmorganlang/bin/SPAdes-3.6.0-Linux/bin/spades.py"]
    # spades_command += ["-1", forward]
    # spades_command += ["-2", reverse]
    # spades_command += ["-s", args.output_dir[sample_id] + "singletons.fastq"]
    # spades_command += ["--careful"]
    # spades_command += ["--memory", str(25)]
    # spades_command += ["--threads", str(args.threads)]
    # spades_command += ["--cov-cutoff", str(10)]
    # spades_command += ["-k", "71,81,91,101,111,121"]
    # spades_command += ["-o", args.output_dir[sample_id] + "assembly"]
    # p_spades = subprocess.Popen(' '.join(spades_command), shell=True, preexec_fn=os.setsid)
    # p_spades.wait()

    megahit_command = [args.executables["megahit"]]
    megahit_command += ["--12", pe]
    megahit_command += ["--read", se]
    megahit_command += ["--k-min", str(k_min), "--k-max", str(k_max)]
    megahit_command += ["--min-count", str(10)]
    megahit_command += ["--k-step", str(10)]
    megahit_command += ["--merge-level", "20,0.98"]
    megahit_command += ["--memory", str(0.25)]
    megahit_output_dir = args.output_dir[sample_id] + "assembly"
    megahit_command += ["--out-dir", megahit_output_dir]
    megahit_command += ["--num-cpu-threads", args.threads]
    megahit_command.append("--verbose")
    megahit_command.append("--no-mercy")  # Recommended for high-coverage datasets
    megahit_command += ["1>", "/dev/null", "2>", "/dev/null"]
    # Run the megahit command using subprocess
    p_megahit = subprocess.Popen(' '.join(megahit_command), shell=True, preexec_fn=os.setsid)
    p_megahit.wait()
    if p_megahit.returncode != 0:
        stdprint("MEGAHIT did not complete successfully!", "err", "\n")

    # Clean-up intermediate output files from MEGAHIT
    stdprint("Cleaning up assembly outputs... ", "out", "")
    contigs_file = args.output_dir[sample_id] + "assembly" + os.sep + "final.contigs.fa"
    if not os.path.isfile(contigs_file):
        sys.exit("ERROR: MEGAHIT assembly was not created! Check " + args.output_dir[sample_id] + "assembly" + os.sep +
                 "log for an error.")
    shutil.move(contigs_file,
                args.output_dir[sample_id] + os.sep + sample_id + "_contigs.fasta")
    shutil.move(args.output_dir[sample_id] + "assembly" + os.sep + "log",
                args.output_dir[sample_id] + os.sep + "megahit.log")
    shutil.rmtree(args.output_dir[sample_id] + "assembly" + os.sep)
    stdprint("done.", "out", "\n")
    return


def prep_reads_for_assembly(sample_id, args, trimmed_reads):
    stdprint("Preparing quality FASTQ files for assembly... ", "out")
    paired_files = list()
    orphan_files = list()
    for fastq in trimmed_reads:
        if re.search(r'pe.[1-2].fq$', fastq):
            paired_files.append(fastq)
        else:
            orphan_files.append(fastq)
    paired_merge = [args.executables["fq2fa"], "--merge"]
    paired_merge += paired_files
    paired_merge.append(args.output_dir[sample_id] + sample_id + "_paired.fa")
    p_pemerge = subprocess.Popen(' '.join(paired_merge), shell=True, preexec_fn=os.setsid)
    p_pemerge.wait()

    orphans = open(args.output_dir[sample_id] + "singletons.fastq", 'w')
    for unpaired_file in orphan_files:
        with open(unpaired_file) as orphan_file:
            for line in orphan_file:
                orphans.write(line)
    orphans.close()
    unpaired_fq2fa = [args.executables["fq2fa"],
                      args.output_dir[sample_id] + "singletons.fastq",
                      args.output_dir[sample_id] + sample_id + "_unpaired.fa"]
    p_sefq2fa = subprocess.Popen(' '.join(unpaired_fq2fa), shell=True, preexec_fn=os.setsid)
    p_sefq2fa.wait()

    reads = glob.glob(args.output_dir[sample_id] + os.sep + sample_id + "_*paired.fa")

    stdprint("done.", "out", "\n")
    return reads


def find_raw_reads(args, sample_id):
    """
    finds the raw read files in the path and stores the file names in a dictionary.
    :param args:
    :param sample_id:
    :return: a dictionary divided by forward and reverse sense
    """
    raw_reads = dict()
    raw_reads["forward"] = ""
    raw_reads["reverse"] = ""
    forward_reads = glob.glob(args.reads + os.sep + "*fastq")
    forward_reads += glob.glob(args.reads + os.sep + "*fq")
    if len(forward_reads) == 0:
        sys.exit("ERROR: Unable to locate fastq files. Must end in either 'fastq' or 'fq'")
    regex_sample = re.compile(sample_id)
    for fastq in forward_reads:
        if regex_sample.search(os.path.basename(fastq)) and raw_reads["forward"] == "":
            raw_reads["forward"] = fastq
        # Ensure there is a single forward fastq file for sample_id
        elif regex_sample.search(os.path.basename(fastq)) and raw_reads["forward"] != "":
            sys.exit("ERROR: More than 2 fastq files match sample ID string in reads directory! "
                     "Please concatenate those files from the same library and try again.")
    if args.reverse:
        reverse_reads = glob.glob(args.reverse + os.sep + "*fastq")
        reverse_reads += glob.glob(args.reverse + os.sep + "*fq")
        for fastq in reverse_reads:
            if re.search(sample_id, fastq) and raw_reads["reverse"] == "":
                raw_reads["reverse"] = fastq
            # Ensure there is a single reverse fastq file for sample_id
            elif regex_sample.search(os.path.basename(fastq)) and raw_reads["reverse"] != "":
                sys.exit("ERROR: More than 2 fastq files match sample ID string in reverse directory! "
                         "Please concatenate those files from the same library and try again.")

    # Check to make sure there are fastq files
    elif len(raw_reads.values()) == 0:
        sys.exit("ERROR: Unable to locate fastq files " + sample_id)
    return raw_reads


class Miffed:
    def __init__(self, fields):
        self.sample = fields[0]
        self.project = fields[1]
        self.selector = fields[2]
        self.vector = fields[3]
        self.screen = fields[4]
        self.selection = fields[5]
        self.num_fosmids = fields[6]
        self.seq_submission_date = fields[7]
        self.glycerol_plate = fields[8]
        self.seq_center = fields[9]
        self.seq_type = fields[10]
        self.read_length = fields[11]
        self.instrument = fields[12]

    def ensure_completeness(self, miffed):
        if self.sample is None:
            stdprint("WARNING: sample name missing in " + miffed + ". Exiting!", "err", "\n")
            sys.exit()
        if self.project is None:
            stdprint("WARNING: " + self.sample + " is not associated with a project!. Skipping this library!",
                     "err", "\n")
            return False
        if self.selector is None:
            stdprint("WARNING: " + self.sample + " is not associated with a human (selector). Skipping this library!",
                     "err", "\n")
            return False
        if self.vector is None:
            stdprint("WARNING: " + self.sample + " is not associated with a vector (e.g. PCC1). Skipping this library!",
                     "err", "\n")
            return False
        if self.screen != "in silico" and self.screen != "functional":
            stdprint("WARNING: " + self.screen + " was provided for screen in sample " + self.sample +
                     ". Is this correct?", "err", "\n")
        return True


def parse_miffed(args):
    """
    Finds the primary project name, submitter, sequencing platform and other information to help with processing
    and storage of the inputs/outputs.
    :param args: parsed command-line arguments from get_options()
    :return:
    """
    miffed = open(args.miffed, 'r')
    line = miffed.readline()
    libraries = dict()
    args.output_dir = dict()
    # excluded contains the libraries that will not be processed for missing MIFFED or non-empty output directories
    excluded = list()
    while line:
        if not line[0] == '#':
            if ',' not in line:
                sys.exit("ERROR: MIFFED input file doesn't seem to be in csv format!")
            fields = line.strip().split(',')
            if len(fields) != 13:
                raise AssertionError(str(len(fields)) + " fields in " + args.miffed + "\n" + line)
            miffed_entry = Miffed(fields)
            complete = miffed_entry.ensure_completeness(args.miffed)
            if complete is False:
                excluded.append(miffed_entry.sample)
            libraries[miffed_entry.sample] = miffed_entry
            args.output_dir[miffed_entry.sample] = args.fabfos_path + os.sep + miffed_entry.project + \
                                                   os.sep + miffed_entry.sample + os.sep
            # Check to ensure the output directory exist or ensure it is empty:
            if os.path.isdir(args.output_dir[miffed_entry.sample]):
                files = glob.glob(args.output_dir[miffed_entry.sample] + "*")
                if len(files) != 0:
                    stdprint("WARNING: output path for " + str(miffed_entry.sample) +
                             " is not empty! Skipping this library.",
                             "err", "\n")
                    excluded.append(miffed_entry.sample)
            else:
                os.makedirs(args.output_dir[miffed_entry.sample])
            project_metadata_file = args.fabfos_path + os.sep + miffed_entry.project + os.sep + \
                                    "FabFos_" + miffed_entry.project + "_metadata.tsv"
            if not os.path.isfile(project_metadata_file):
                project_metadata = open(project_metadata_file, 'w')
                project_metadata.write(args.metadata_header)
                project_metadata.close()
        line = miffed.readline()
    miffed.close()
    for library in excluded:
        # if library != "Zach-s-Beaver-positives-plate-1_S2":
        #     args.output_dir.pop(library)
        args.output_dir.pop(library)
    return libraries, args


def clean_intermediates(sample_id, args):
    """
    Function removes largely useless alignment files
    :param sample_id: The unique identifier for the sample - used for output directory creation
    :param args: parsed command-line arguments from get_options()
    :return:
    """
    # Remove the alignment files
    working_dir = args.output_dir[sample_id] + os.sep
    os.remove(working_dir + sample_id + ".sam")
    os.remove(working_dir + sample_id + "_sorted.bam")

    # Remove bam2fastq outputs:
    for std_file in glob.glob(working_dir + "bam2fastq*"):
        os.remove(std_file)

    # Remove trimmomatic outputs
    trimmomatic_outputs = glob.glob(working_dir + sample_id + "_trim*fq")
    for trim_out in trimmomatic_outputs:
        os.remove(trim_out)

    os.remove(working_dir + "singletons.fastq")
    return


# Heng Li's readfq generator function for reading fastq files
def readfq(fp):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, length, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                length += len(l) - 1
                if length >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


def find_num_reads(file_list):
    num_reads = 0
    for fastq in file_list:
        fq_in = open(fastq, 'r')
        for name, seq, qual in readfq(fq_in):
            num_reads += 1
        fq_in.close()
    return num_reads


def map_ends(args, sample_id):
    """
    Function for using blastn to align fosmid ends in a file to find the well of each contig
    """
    # Index the contigs
    sample_prefix = args.output_dir[sample_id] + os.sep + sample_id
    blastdb_command = [args.executables["makeblastdb"]]
    blastdb_command += ["-in", sample_prefix + "_contigs.fasta"]
    blastdb_command += ["-dbtype", "nucl"]
    blastdb_command += ["-out", sample_prefix]
    blastdb_command += ["1>/dev/null", "2>/dev/null"]
    p_makeblastdb = subprocess.Popen(" ".join(blastdb_command), shell=True, preexec_fn=os.setsid)
    p_makeblastdb.wait()

    stdprint("Aligning fosmid ends to assembly... ", "out")
    blastn_command = [args.executables["blastn"]]
    blastn_command += ["-db", sample_prefix]
    blastn_command += ["-outfmt", "\"6",
                       "sseqid", "slen", "qseqid", "pident", "length", "sstrand", "sstart", "send", "bitscore"
                       "\""]
    blastn_command += ["-query", args.ends]
    blastn_command += ["-perc_identity", str(95)]
    blastn_command += ["-out", sample_prefix + "_endsMapped.tbl"]
    blastn_command += ["1>/dev/null", "2>/dev/null"]
    p_align_ends = subprocess.Popen(" ".join(blastn_command), shell=True, preexec_fn=os.setsid)
    p_align_ends.wait()
    stdprint("done", "out", "\n")

    return


def get_fosmid_end_name(string):
    full_name = string.split('.')[1]
    prefix = string.split('.')[0]
    # print "string= ", string, "\tfull= ", full_name, "\tprefix= ", prefix
    direction, name = full_name.split('_')
    name = prefix + '.' + name
    return name


class EndAlignment:
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
        self.name = ""
        self.direction = ""

    def parse_fosmid_end_name(self):
        full_name = self.qseqid.split('.')[1]
        prefix = self.qseqid.split('.')[0]
        direction, name = full_name.split('_')
        self.direction = direction[-1]
        self.name = prefix + '.' + name

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


def parse_end_alignments(args, sample_id, fosmid_assembly, ends_stats):
    """
    Function parses BLAST alignments of fosmid ends mapped to the MEGAHIT assemblies
    :param args: parsed command-line arguments from get_options()
    :param sample_id: string representing the current sample name being processed
    :param fosmid_assembly: dictionary with headers as keys and sequences for values
    :param ends_stats: dictionary containing information on the fosmid ends
    :return: dictionary of contig names (keys) and list of EndAlignment objects (values) or empty list
    """
    ends_mapping = dict()
    ends_stats["aligned"] = set()
    ends_stats["unaligned"] = set()
    blast_output = args.output_dir[sample_id] + os.sep + sample_id + "_endsMapped.tbl"
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
        if alignment.length > 200 and alignment.bitscore > 500:
            alignment.parse_fosmid_end_name()
            if alignment.name not in ends_stats["failed"]:
                ends_stats["aligned"].add(alignment.name)

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

    unaligned = ends_stats["all_clones"].difference(ends_stats["aligned"])
    for fosmid_end in unaligned:
        if fosmid_end not in ends_stats["failed"]:
            ends_stats["unaligned"].add(fosmid_end)
    ends_stats["Num_unaligned"] = len(ends_stats["unaligned"])
    ends_stats["Num_aligned"] = len(ends_stats["aligned"])

    return ends_mapping, ends_stats


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

    if not forward_aligned:
        stdprint("WARNING: Forward strand did not align in single-pair clone " + clone, "err", "\n")
    if not reverse_aligned:
        stdprint("WARNING: Reverse strand did not align in single-pair clone " + clone, "err", "\n")
    if len(positions_dict["minus_strand"]) == 0 or len(positions_dict["plus_strand"]) == 0:
        stdprint("WARNING: Forward and reverse mate-pairs of " + clone + " aligned to the same strand of " + contig +
                 "!\nThis sequence should not be trusted due to assembly or library preparation errors", "err", "\n")

    positions = positions_dict["plus_strand"] + positions_dict["minus_strand"]
    start = min(positions)
    end = max(positions)
    return start, end


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


def assign_clones(ends_mapping, ends_stats, fosmid_assembly):
    """
    Assigns fosmid-end sequences to assembled contigs in fosmid_assembly
    :param ends_mapping:
    :param ends_stats:
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
    ends_stats["Single-orphan"] = len(single_orphans)
    ends_stats["Unassigned"] = len(unassigned_contigs)
    ends_stats["Multi-pair"] = len(multi_pairs)
    ends_stats["Single-pair"] = len(single_pairs)
    ends_stats["Multi-orphan"] = len(multi_orphans)
    return clone_map, ends_stats, multi_fosmid_map


def prune_and_scaffold_fosmids(args, sample_id, clone_map, multi_fosmid_map):
    fragments_tsv = args.output_dir[sample_id] + os.sep + "potential_fosmid_fragments.tsv"

    try:
        fragments = open(fragments_tsv, 'w')
    except IOError:
        sys.exit("Unable to open " + fragments_tsv + " for writing!")

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
                            new_fosmid = FosmidClone(leading_node + ',' + trailing_node, prefix + suffix)
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
            # Not confident in scaffolding these contigs
            # Write the FosmidClone information to a potential_fosmid_fragments.tsv file
            for fosmid in multi_fosmid_map[clone]:
                fragments.write(fosmid.get_info() + "\n")

    fragments.close()

    return clone_map


def write_fosmid_assignments(args, sample_id, clone_map):
    coord_name = args.output_dir[sample_id] + os.sep + sample_id + "_fosmid_coordinates.tsv"
    fasta_name = args.output_dir[sample_id] + os.sep + sample_id + "_formatted.fasta"
    try:
        fos_coord = open(coord_name, 'w')
    except IOError:
        sys.exit("Unable to open " + coord_name + " for writing!")

    try:
        fos_fasta = open(fasta_name, 'w')
    except IOError:
        sys.exit("Unable to open " + fasta_name + " for writing!")

    fos_coord.write("#Fosmid\tContig\tEvidence\tEnds\tSequence\n")
    for fosmid_clone in clone_map:
        fos_coord.write(fosmid_clone.clone + "\t" +
                        fosmid_clone.contig + "\t" +
                        fosmid_clone.evidence + "\t" +
                        fosmid_clone.ends + "\t" +
                        fosmid_clone.sequence + "\n")
        fos_fasta.write(">" + fosmid_clone.clone + " " + fosmid_clone.contig[1:] + " " + fosmid_clone.evidence + "\n")
        fos_fasta.write(fosmid_clone.sequence + "\n")

    fos_coord.close()
    fos_fasta.close()

    return


def read_fasta(fasta):
    """
    Function to read a fasta file and store it in a dictionary. Headers are keys, sequences are values
    :param fasta: A fasta file name (and path, if necessary)
    :return: dictionary with headers as keys and sequences for values
    """
    fasta_dict = dict()
    header = ""
    sequence = ""
    with open(fasta) as assembly:
        line = assembly.readline()
        while line:
            line = line.strip()
            if line[0] == '>':
                if header != "":
                    fasta_dict[header] = sequence
                # Takes the first word in the header (to match alignment outputs)
                header = line.split(" ")[0]
                sequence = ""
            else:
                sequence += line
            line = assembly.readline()
        fasta_dict[header] = sequence

    return fasta_dict


def get_assembly_stats(sample_id, args, fosmid_assembly):
    assembly_stats = dict()
    sample_prefix = args.output_dir[sample_id] + os.sep + sample_id
    nx_command = [args.executables["getNx"]]
    nx_command += ["-i", sample_prefix + "_contigs.fasta"]
    nx_command += ["-o", sample_prefix + "_nx.csv"]
    nx_command += ["-m", str(2000)]
    nx_command.append(">" + args.output_dir[sample_id] + os.sep + sample_id + "_nx.txt")

    p_nx = subprocess.Popen(' '.join(nx_command), shell=True, preexec_fn=os.setsid)
    p_nx.wait()
    if p_nx.returncode != 0:
        stdprint("getNx did not finish successfully!", "err", "\n")
        sys.exit()

    try:
        nx_stats = open(sample_prefix + "_nx.csv", 'r')
    except IOError:
        sys.exit("getNx did not finish successfully: unable to open " + sample_prefix + "_nx.csv for reading!")

    nx_stats.readline()
    line = nx_stats.readline()
    while line:
        line = line.strip()
        proportion, n = line.split(',')
        if float(proportion) == 0.5:
            stdprint("N50 = " + str(n), "out", "\n")
            assembly_stats["N50"] = str(n)
        if float(proportion) == 0.9:
            assembly_stats["N90"] = str(n)
        if float(proportion) == 0:
            stdprint("Longest contig = " + str(n) + "bp", "out", "\t")
        line = nx_stats.readline()

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

    os.remove(args.output_dir[sample_id] + os.sep + sample_id + "_nx.txt")

    return assembly_stats


def update_metadata(metadata_file, sample_id, library_metadata, read_stats, assembly_stats, ends_stats):
    try:
        metadata = open(metadata_file, 'a')
    except IOError:
        sys.exit("ERROR: Unable to open " + metadata_file + " to append library metadata!")
    metadata.write(sample_id + "\t")
    metadata.write(library_metadata.project + "\t")
    metadata.write(library_metadata.selector + "\t")
    metadata.write(library_metadata.vector + "\t")
    metadata.write(library_metadata.screen + "\t")
    metadata.write(library_metadata.selection + "\t")
    metadata.write(library_metadata.num_fosmids + "\t")
    metadata.write(library_metadata.seq_submission_date + "\t")
    metadata.write(library_metadata.glycerol_plate + "\t")
    metadata.write(library_metadata.seq_center + "\t")
    metadata.write(library_metadata.seq_type + "\t")
    metadata.write(library_metadata.read_length + "\t")
    metadata.write(library_metadata.instrument + "\t")
    metadata.write(strftime("%Y-%m-%d") + "\t")
    metadata.write(read_stats["num_raw_reads"] + "\t")
    metadata.write(read_stats["percent_off_target"] + "\t")
    metadata.write(str(read_stats["number_trimmed_reads"]) + "\t")
    metadata.write(assembly_stats["Contigs"] + "\t")
    metadata.write(assembly_stats["N50"] + "\t")
    metadata.write(assembly_stats["N90"] + "\t")
    metadata.write(assembly_stats["27kbp"] + "\t")
    metadata.write(assembly_stats["50kbp"] + "\t")
    metadata.write(str(ends_stats["Single-pair"]) + "\t")
    metadata.write(str(ends_stats["Multi-pair"]) + "\t")
    metadata.write(str(ends_stats["Single-orphan"]) + "\t")
    metadata.write(str(ends_stats["Multi-orphan"]) + "\t")
    metadata.write(str(ends_stats["Unassigned"]) + "\t")
    metadata.write(str(ends_stats["Num_total"]) + "\t")
    metadata.write(str(ends_stats["Num_aligned"]) + "\t")
    metadata.write(str(ends_stats["Num_unaligned"]) + "\t")
    metadata.write(str(ends_stats["Num_failed"]) + "\n")

    metadata.close()
    return


def get_fosmid_ends_stats(ends_fasta):
    ends_stats = dict()
    ends_stats["failed"] = set()
    ends_stats["all_clones"] = set()
    for clone in ends_fasta:
        clone_name = get_fosmid_end_name(clone[1:])
        ends_stats["all_clones"].add(clone_name)
        if len(ends_fasta[clone]) < 100:
            ends_stats["failed"].add(clone_name)

    ends_stats["Num_total"] = len(ends_stats["all_clones"])
    ends_stats["Num_failed"] = len(ends_stats["failed"])
    return ends_stats


def write_fosmid_end_failures(args, sample_id, ends_stats):
    """
    Writes the names of fosmid ends that could not be aligned and failed
    :param args:
    :param sample_id:
    :param ends_stats:
    :return:
    """
    sample_prefix = args.output_dir[sample_id] + os.sep + sample_id
    fosmid_end_failures = sample_prefix + "_fosmid_end_failures.tsv"

    try:
        failure_file = open(fosmid_end_failures, 'w')
    except:
        raise IOError("ERROR: cannot open file: " + fosmid_end_failures)
    failure_file.write("#Fosmid-end\tcategory\n")
    for end in ends_stats["failed"]:
        failure_file.write(end + "\tSequencing failure\n")
    for end in ends_stats["unaligned"]:
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


def main():
    # TODO: Include an init function to initialize a new directory to be a master FabFos repository
    args = get_options()
    args = review_arguments(args)
    args = find_executables(args)
    libraries, args = parse_miffed(args)

    ends_stats = dict()
    if args.ends:
        ends_fasta = read_fasta(args.ends)
        ends_stats = get_fosmid_ends_stats(ends_fasta)
    else:
        ends_stats["Single-pair"] = 0
        ends_stats["Multi-pair"] = 0
        ends_stats["Single-orphan"] = 0
        ends_stats["Multi-orphan"] = 0
        ends_stats["Unassigned"] = 0
        ends_stats["Total_ends"] = 0
        ends_stats["Num_missing"] = 0
        ends_stats["Num_stunted"] = 0

    # Sequence processing and filtering:
    for sample_id in args.output_dir.keys():
        read_stats = dict()
        stdprint("Processing raw data for " + sample_id, "out", "\n")
        stdprint("Outputs for " + sample_id + " will be found in " + args.output_dir[sample_id], "out", "\n")
        raw_reads = find_raw_reads(args, sample_id)
        num_raw_reads = find_num_reads(raw_reads.values())
        read_stats["num_raw_reads"] = str(num_raw_reads)
        stdprint("Number of raw reads = " + str(num_raw_reads), "out", "\n")
        filtered_reads = filter_backbone(sample_id, args, raw_reads)
        num_filtered_reads = find_num_reads(filtered_reads)
        stdprint(str(num_raw_reads - num_filtered_reads) + " reads removed from backbone filtering (" +
                 str(((num_raw_reads - num_filtered_reads)*100)/num_raw_reads) + "%).",
                 "out",
                 "\n")
        read_stats["percent_off_target"] = str(((num_raw_reads - num_filtered_reads)*100)/num_raw_reads)
        trimmed_reads = quality_trimming(sample_id, args, filtered_reads)
        num_reads_assembly = find_num_reads(trimmed_reads)
        read_stats["number_trimmed_reads"] = str((num_filtered_reads - num_reads_assembly)*100/num_reads_assembly)
        reads = prep_reads_for_assembly(sample_id, args, trimmed_reads)
        # TODO: Include an optional minimus2 module to further assemble the contigs
        assemble_fosmids(sample_id, args, reads, 71, 241, trimmed_reads)
        fosmid_fasta = args.output_dir[sample_id] + os.sep + sample_id + "_contigs.fasta"
        fosmid_assembly = read_fasta(fosmid_fasta)
        assembly_stats = get_assembly_stats(sample_id, args, fosmid_assembly)
        # For mapping fosmid ends:
        if args.ends:
            map_ends(args, sample_id)
            ends_mapping, ends_stats = parse_end_alignments(args, sample_id, fosmid_assembly, ends_stats)
            clone_map, ends_stats, multi_fosmid_map = assign_clones(ends_mapping, ends_stats, fosmid_assembly)
            clone_map = prune_and_scaffold_fosmids(args, sample_id, clone_map, multi_fosmid_map)
            write_fosmid_assignments(args, sample_id, clone_map)
            write_fosmid_end_failures(args, sample_id, ends_stats)
        clean_intermediates(sample_id, args)
        project_metadata_file = args.fabfos_path + os.sep + libraries[sample_id].project + os.sep + \
                                "FabFos_" + libraries[sample_id].project + "_metadata.tsv"
        update_metadata(project_metadata_file, sample_id, libraries[sample_id], read_stats, assembly_stats, ends_stats)
        update_metadata(args.master_metadata, sample_id, libraries[sample_id], read_stats, assembly_stats, ends_stats)


main()
