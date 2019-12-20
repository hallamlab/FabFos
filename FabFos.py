#!/usr/bin/env python3

try:
    import argparse
    import sys
    import glob
    import os
    import traceback
    import re
    import subprocess
    import shutil
    import logging
    from time import strftime
except ImportWarning:
    sys.stderr.write("Could not load some user defined module functions")
    sys.stderr.write(traceback.print_exc(10))
    sys.exit(3)

"""
FabFos: a pipeline for automatically performing quality controls, assembly, and data storage management
for fosmid sequence information. Circa 2016 - Hallam Lab, UBC
"""

__version__ = "1.1"
__author__ = "Connor Morgan-Lang"
__license__ = "GPL-v3"
__maintainer__ = "Connor Morgan-Lang"
__email__ = "c.morganlang@gmail.com"
__status__ = "Unstable"

#################################### Classes begin ################################################

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


class Sample:
    def __init__(self, sample_id):
        # General sample information
        self.id = sample_id
        self.output_dir = ""
        self.assembled_fosmids = ""
        self.forward_reads = None
        self.reverse_reads = None

        # Control flow
        self.overwrite = False
        self.assemble = True
        self.map_ends = True
        self.exclude = False

        # Nanopore-specific
        self.nanopore = ""
        self.error_correction = True


class Miffed(Sample):
    def __init__(self, sample_id):
        Sample.__init__(self, sample_id)
        self.project = ""
        self.selector = ""
        self.vector = ""
        self.screen = ""
        self.selection = ""
        self.num_fosmids = ""
        self.seq_submission_date = ""
        self.glycerol_plate = ""
        self.seq_center = ""
        self.seq_type = ""
        self.read_length = ""
        self.instrument = ""

    def populate_info(self, fields):
        self.project = fields[1]
        self.selector = fields[2]
        self.vector = fields[3]
        self.screen = fields[4]
        self.selection = fields[5]
        self.seq_submission_date = fields[7]
        self.glycerol_plate = fields[8]
        self.seq_center = fields[9]
        self.seq_type = fields[10]
        self.instrument = fields[12]
        try:
            self.num_fosmids = int(fields[6])
        except ValueError:
            logging.error("Number of fosmids field (column 7) is not an integer!\n")
            sys.exit(3)
        try:
            self.read_length = int(fields[11])
        except ValueError:
            logging.error("Read length field (column 12) is not an integer!\n")
            sys.exit(3)

    def ensure_completeness(self, miffed):
        if self.id is None:
            logging.warning("Sample name missing in " + miffed + ". Exiting!\n")
            sys.exit()
        if self.project is None:
            logging.warning(self.id + " is not associated with a project!. Skipping this library!\n")
            return False
        if self.selector is None:
            logging.warning(self.id + " is not associated with a human (selector). Skipping this library!\n")
            return False
        if self.vector is None:
            logging.warning(self.id + " is not associated with a vector (e.g. PCC1). Skipping this library!\n")
            return False
        if self.screen != "in silico" and self.screen != "functional":
            logging.warning(self.screen + " was provided for screen in sample " + self.id + ". Is this correct?\n")
        return True


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
        # full_name = self.qseqid.split('.')[1]
        prefix = self.qseqid.split('.')[0]
        direction, name = prefix.split('_')
        self.direction = direction[-1]
        self.name = prefix + '.' + name


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

#################################### Classes end ################################################

def remove_singletons(clusters):
    for c in clusters:
        if "size=1;" not in c.id:
            cluster = re.findall('\w+', c.id)[0]
            size = re.findall('(?<=size\=)\w+', c.id)[0]
            output_row = [cluster, size, c.seq, c.seq.reverse_complement()]
            final_output.append(output_row)
    return final_output

def bam2fastq(input_file, sample_id, output_dir, executables):
    fastq_extract = [executables["bam2fastq"], "-o",
                     output_dir + os.sep + sample_id + input_file + ".fasta", "--no-aligned", input_file]
    fastq_extract += ["1>", output_dir + os.sep + "bam2fastq.stdout"]
    fastq_extract += ["2>", output_dir + os.sep + "bam2fastq.stderr"]
    p_bam2fastq = subprocess.Popen(' '.join(fastq_extract), shell=True, preexec_fn=os.setsid)
    p_bam2fastq.wait()
    filtered_reads = glob.glob(output_dir + os.sep + sample_id + input_file + "*fastq")
    logging.info("done.\n")
    return filtered_reads


def get_dict_keys(dictionary):
    """
    A function for returning the sorted list of dictionary keys, regardless of the Python version.
    :param dictionary: dict() object
    :return: list() of sorted dictionary keys
    """
    dictionary.keys()
    return


def get_options():
    parser = argparse.ArgumentParser(description="Pipeline for processing and organizing fosmid sequence information.\n"
                                                 "NOTE: a maximum of 2 sequence files are permitted.",add_help=False)
    reqs = parser.add_argument_group(title="Required arguments")
    reqs.add_argument("-m", "--miffed", type=str, required=True,
                      help="The minimum information for fosmid environmental data (e.g., sample ID, "
                           "sequencing platform, environment, project) in a comma-separated file. [.csv]")
    reqs.add_argument("-b", "--background", type=str, required=True,
                      help="Path to the fosmid background database [.fasta]")
    reqs.add_argument("-r", "--reads", type=str, required=True,
                      help="Path to the sequenced reads directory. Can be in FastQ or BAM format. This parameter is for the
                           "single-end reads Sanger), the forward strand file, or the interleaved paired-end file. [.fastq]")

    nanopore = parser.add_argument_group(title="Nanopore-specific [development] options")
    nanopore.add_argument("--nanopore_reads", help="A FASTA file containing nanopore reads to be used in assembly.",
                          type=str, required=False)
    nanopore.add_argument("--skip_correction", help="Do not perform error-correction of nanopore reads using proovread",
                          action="store_true", default=False, required=False)

    opts = parser.add_argument_group(title="Optional arguments")

    opts.add_argument("-c", "--bb_length", required=False, default= str(40),
                      help="Enter the minimum length of the clustered reads that must align to the edge of the backbone")
    opts.add_argument("-l", "--cluster_length", required=False, default= str(100),
                      help="Enter the length of the trimmed sequences you want clustered.")
    opts.add_argument("-n", "--backbone", required=False, default="pCC1fos",
                      help="Enter the name of the vector backbone you're using")
    opts.add_argument("-t", "--type", choices=["B", "F"], required=False, default="F",
                      help="Enter B if input type is BAM, F for FastQ")
    opts.add_argument("-2", "--reverse", type=str, required=False,
                      help="Path to the directory containing reverse-end reads (if applicable) [.fastq]")
    opts.add_argument("-i", "--interleaved", required=False, default=False, action="store_true",
                      help="Flag indicating the reads are interleaved "
                           "(i.e. forward and reverse pairs are in the same file)")
    opts.add_argument("-a", "--assembler", choices=["spades", "megahit"], required=False, default="spades",
                      help="Genome assembly software to use.")
    opts.add_argument("-f", "--fabfos_path", type=str, required=False,
                      default="/mnt/nfs/sharknado/LimsData/FabFos/",
                      help="Path to FabFos database on sharknado [DEFAULT = /mnt/nfs/sharknado/LimsData/FabFos/]")
    opts.add_argument("-T", "--threads", type=str, required=False, default=str(8),
                      help="The number of threads that can be used [DEFAULT = 8]")
    opts.add_argument("-e", "--ends", required=False, default=None,
                      help="FASTA file containing fosmid ends - these will be used for alignment to fosmid contigs."
                           " [.fasta]")
    opts.add_argument("-v", "--verbose", required=False, default=False, action="store_true",
                      help="Increase the level of verbosity in runtime log.")
    opts.add_argument("-h", "--help",
                      action="help", help="Show this help message and exit")
    
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


def prep_logging(log_file_name, verbosity):
    logging.basicConfig(level=logging.DEBUG,
                        filename=log_file_name,
                        filemode='w',
                        datefmt="%d/%m %H:%M:%S",
                        format="%(asctime)s %(levelname)s:\n%(message)s")
    if verbosity:
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


def review_arguments(args):
    """
    Function that ensures the proper information has been provided and adds system information to args
    :param args: parsed command-line arguments from get_options()
    :return: args with more information
    """

    # Setup the global logger and main log file
    log_file_name = args.fabfos_path + os.sep + "FabFos_log.txt"
    prep_logging(log_file_name, args.verbose)

    # Review the provided arguments:
    if not os.path.isdir(args.reads):
        logging.error(args.reads + " is not a valid directory!\n")
        sys.exit(3)
    if args.reverse and not os.path.isdir(args.reverse):
        logging.error(args.reverse + " is not a valid directory!\n")
        sys.exit(3)

    if not os.path.isfile(args.background):
        logging.error(args.background + " does not exist!\n")
        sys.exit(3)

    if args.reverse and args.interleaved:
        logging.error("Reads cannot be interleaved and also have separate forward- and reverse-FASTQ files!\n")
        sys.exit(3)

    # Add new information:
    args.os = os_type()
    if sys.version_info > (2, 9):
        args.py_version = 3
    else:
        args.py_version = 2

    args.adaptor_trim = "/usr/local/share/Trimmomatic-0.35/adapters/"

    args.metadata_header = "#Sample Name (LLLLL-PP-WWW)\tProject\tHuman selector\tVector Name\t" \
                           "Screen [in silico | functional]\tSelection criteria\tNumber of fosmids\t" \
                           "Sequencing submission date (YYYY-MM-DD)\tGlycerol plate name\tSequencing center\t" \
                           "Sequencing type\tRead length\tInstrument\tDate of FabFos analysis\tNumber of reads\t" \
                           "% off-target reads\tNumber of trimmed reads\tAssembler version\tNumber of Contigs\t" \
                           "N50\tN90\tContigs > 27kbp\tContigs > 50kbp\t" \
                           "# single-pair\t# multi-pair\t# single-orphan\t# multi-orphan\t# Unidentifiable\t" \
                           "# Fosmid ends(X2 = #reads)\t# aligned fosmid ends\t# unaligned fosmid ends\t# Failed ends\n"

    if not os.path.isfile(args.fabfos_path + os.sep + "FabFos_master_metadata.tsv"):
        logging.error("This fabfos_path directory does not contain a FabFos_master_metadata.tsv file. "
                      "Are you sure its a bona fide FabFos repository?")
        sys.exit(3)

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
    # TODO: Write the tool versions to the log file
    required_execs = ["bwa", "samtools", "megahit", "bam2fastq", "trimmomatic", "fq2fa", "blastn", "makeblastdb",
                      "splitFASTA", "getNx", "proovread", "canu"]
    args.executables = dict()
    for executable in required_execs:
        args.executables[executable] = which(executable)
        if args.executables[executable] is None:
            raise EnvironmentError("Unable to find executable for " + executable)
    return args


def check_index(path, bwa_path):
    extensions = ['.bwt', '.pac', '.ann', '.amb', '.sa']
    index_parts = [path+e for e in extensions]
    for part in index_parts:
        if not os.path.isfile(part):
            logging.info("\nUnable to find BWA index for " + path + ". Indexing now... ")
            index_command = [bwa_path, "index", path]
            index_dir = os.path.dirname(os.path.abspath(path))
            index_command += ["1>", index_dir + os.sep + "bwa_index.stdout"]
            index_command += ["2>", index_dir + os.sep + "bwa_index.stderr"]
            p_index = subprocess.Popen(' '.join(index_command), shell=True, preexec_fn=os.setsid)
            p_index.wait()
            if p_index.returncode != 0:
                sys.exit("ERROR: bwa index did not complete successfully")
            logging.info("done.\nResuming alignment... ")
            break
    return


def filter_backbone(sample, args, raw_reads):
    """
    Function to generate fastq files that do not contain any sequences in `background`.
    Depends on bwa, samtools, and bam2fastq
    :param sample: Miffed object with information of current sample
    :param args: parsed command-line arguments from get_options()
    :param raw_reads: Dictionary containing forward and reverse fastq file names
    :return: list of fastq files containing the filtered reads
    """
    logging.info("Filtering off-target reads... ")
    check_index(args.background, args.executables["bwa"])
    align_command = [args.executables["bwa"], "mem", "-t", args.threads, args.background, raw_reads["forward"]]
    if args.reverse or args.interleaved:
        align_command.append(raw_reads["reverse"])
    align_command.append("1>")
    sam_file = sample.output_dir + sample.id + ".sam"
    align_command.append(sam_file)
    align_command += ["2>", "/dev/null"]
    p_mem = subprocess.Popen(' '.join(align_command), shell=True, preexec_fn=os.setsid)
    p_mem.wait()

    # Use samtools to convert sam to bam
    bam_convert = [args.executables["samtools"], "view", "-bS", "-@", args.threads, sam_file, "|"]
    # Piping to samtools sort
    bam_file = sample.output_dir + os.sep + sample.id + "_sorted.bam"
    bam_convert.append(args.executables["samtools"])
    bam_convert += ["sort", "-@", args.threads, "-"]
    bam_convert += ["-o", bam_file]
    p_samtools_stdout = open(sample.output_dir + os.sep + "samtools.stdout", 'w')
    p_samtools_stderr = open(sample.output_dir + os.sep + "samtools.stderr", 'w')
    p_samtools = subprocess.Popen(' '.join(bam_convert), shell=True, preexec_fn=os.setsid,
                                  stdout=p_samtools_stdout, stderr=p_samtools_stderr)
    p_samtools.wait()
    if p_samtools.returncode != 0:
        logging.error("BAM file was not successfully created and sorted!\n")
        sys.exit(3)
    p_samtools_stdout.close()
    p_samtools_stderr.close()

    # extract the unaligned reads using bam2fastq
    filtered_reads = bam2fastq(bam_file, sample.id, sample.output_dir, args.executables)
    return filtered_reads, bam_file

#bam_file: bam file containing the alignment map 
#output_dir: Directory where all the files are to be output
#bl: Length of the backbone to which the reads are to be aligned
#rl: Length of the reads that need to be clustered

#bam_file: bam file containing the alignment map 
#output_dir: Directory where all the files are to be output
#bl: minimum length of the backbone to which the reads are to be aligned
#rl: Length of the reads that need to be clustered

def clustering(bam_file, output_dir, bl, rl):

    backbone_name = args.backbone

    command_list = []
    subprocess.call("mkdir clustering_log", shell=True).

    #Isolate the reads that map to the first and last bl bp of the backbone
    sam_view1 = ['samtools', 'view', '-F 4',
                     str(bam_file),
                     str(backbone_name +':0' +'-'+ str(bl)),
                     '-o',
                     str(clustering_log + os.sep + '5_040.bam')]
    sam_view2 = ['samtools', 'view', '-F 4',
                     str(bam_file),
                     str(backbone_name + ':' + str(len(backbone_fasta.seq)-bl) +'-'+ str(len(backbone_fasta.seq))),
                     '-o',
                     str(clustering_log + os.sep + '3_040.bam')]

    command_list.append(sam_view1)
    command_list.append(sam_view2)

    #Isolate the fwd reads that map to the first and last bl bp of the backbone and create a new file for them in output directory
    sam_view3 = ['samtools', 'view', '-F', '0x10',
                     str(clustering_log + os.sep + '5_040.bam'),
                     '-o',
                     str(clustering_log + os.sep +'040-fwd.bam')]
    sam_view4 = ['samtools', 'view', '-F', '0x10',
                     str(clustering_log + os.sep + '3_040.bam'),
                     '-o',
                     str(clustering_log + os.sep +'l40-fwd.bam')]

    command_list.append(sam_view3)
    command_list.append(sam_view4)

    #Isolate the rev reads that map to the first and last bl bp of the backbone and create a new file for them in output directory
    sam_view5 = ['samtools', 'view', '-f', '0x10',
                     str(clustering_log + os.sep + '5_040.bam'),
                     '-o',
                     str(clustering_log + os.sep +'040-rev.bam')]
    sam_view6 = ['samtools', 'view', '-f', '0x10',
                     str(clustering_log + os.sep + '3_040.bam'),
                     '-o',
                     str(clustering_log + os.sep +'l40-rev.bam')]

    command_list.append(sam_view5)
    command_list.append(sam_view6)

    # Generate .fasta file with fwd reads
    sam_fasta1 = ['samtools', 'fasta',
                      str(clustering_log + os.sep +'040-fwd.bam'),
                      '>',
                      str(clustering_log + os.sep +'040-fwd.fasta')]
    sam_fasta2 = ['samtools', 'fasta',
                      str(clustering_log + os.sep +'l40-fwd.bam'),
                      '>',
                      str(clustering_log + os.sep +'l40-fwd.fasta')]

    command_list.append(sam_fasta1)
    command_list.append(sam_fasta2)

    # Generate .fasta file with rev reads
    sam_fasta3 = ['samtools', 'fasta',
                      str(clustering_log + os.sep +'040-rev.bam'),
                      '>',
                      str(clustering_log + os.sep +'040-rev.fasta')]
    sam_fasta4 = ['samtools', 'fasta',
                      str(clustering_log + os.sep +'l40-rev.bam'),
                      '>',
                      str(clustering_log + os.sep +'l40-rev.fasta')]

    command_list.append(sam_fasta3)
    command_list.append(sam_fasta4)

    for command in command_list:
        subprocess.call(' '.join(command), shell = True)

    # Open fwd and rev .fasta files
    fwd1 = SeqIO.parse(str(clustering_log + os.sep+'040-fwd.fasta'), "fasta")

    rev1 = SeqIO.parse(str(clustering_log + os.sep+'040-rev.fasta'), "fasta")

    fwd2 = SeqIO.parse(str(clustering_log + os.sep+'l40-fwd.fasta'), "fasta")

    rev2 = SeqIO.parse(str(clustering_log + os.sep+'l40-rev.fasta'), "fasta")

    # Append fwd reads and reverse_complement of rev reads to recs = []
    recs1 = []
    for f in fwd1:
        recs1.append(f)

    for r in rev1:
        r.seq = r.seq.reverse_complement()
        recs1.append(r)

    recs2 = []
    for f in fwd2:
        recs2.append(f)

    for r in rev2:
        r.seq = r.seq.reverse_complement()
        recs2.append(r)

    # Trim sequences to rl bp from defined site

    final_recs1 = []
    final_recs2 = []

    first_eight = str(backbone_fasta.seq[0:8])
    last_eight = str(backbone_fasta.seq[len(backbone_fasta.seq)-8:len(backbone_fasta.seq)])

    for rec in recs1:
        if first_eight in rec.seq:
            new_seq = rec.seq[0:str(rec.seq).index(first_eight)]
            if len(new_seq) >= rl:
                trim_seq = new_seq[len(new_seq)-rl:len(new_seq)]
                rec.seq = trim_seq
                final_recs1.append(rec)

    for rec in recs2:
        if last_eight in rec.seq:
            new_seq = rec.seq[str(rec.seq).index(last_eight)+8:]
            if len(new_seq) >= rl:
                trim_seq = new_seq[0:rl]
                rec.seq = trim_seq
                final_recs2.append(rec)

    SeqIO.write(final_recs1, str(clustering_log + os.sep +'5-040.fasta'), "fasta")
    SeqIO.write(final_recs2, str(clustering_log + os.sep +'3-040.fasta'), "fasta")

# Cluster trimmed sequences in .fasta
    uclust_cmd1 = ["usearch",
                  "-cluster_fast", clustering_log + os.sep + '5-040.fasta',
                  "-uc", clustering_log + os.sep +'5-040.uc',
                  "-id", str(0.9),
                  "-consout", clustering_log + os.sep +'5-040_cluster.fasta',
                  "-sizeout",
                  "-centroids", clustering_log + os.sep + '5-040_centroids.fasta']

    uclust_cmd2 = ["usearch",
                  "-cluster_fast", clustering_log + os.sep + '3-040.fasta',
                  "-uc", clustering_log + os.sep +'3-040.uc',
                  "-id", str(0.9),
                  "-consout", clustering_log + os.sep +'3-040_cluster.fasta',
                  "-sizeout",
                  "-centroids", clustering_log + os.sep + '3-040_centroids.fasta']

    subprocess.call(' '.join(uclust_cmd1), shell = True)
    subprocess.call(' '.join(uclust_cmd2), shell = True)

# Remove singletons

    clusters1 = SeqIO.parse(str(clustering_log + os.sep +'5-040_cluster.fasta'), "fasta")
    clusters2 = SeqIO.parse(str(clustering_log + os.sep +'3-040_cluster.fasta'), "fasta")

    final_output1 = remove_singletons(clusters1)
    final_output2 = remove_singletons(clusters2)

    with open(ouput_dir + os.sep + 'cluster1.csv', 'w') as writeFile:
        writer = csv.writer(writeFile)
        for fo in final_output1:
            writer.writerow(fo)

    with open(ouput_dir + os.sep + 'cluster2.csv', 'w') as writeFile:
        writer = csv.writer(writeFile)
        for fo in final_output2:
            writer.writerow(fo)


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


def write_new_fasta(fasta_dict, fasta_name, headers=None):
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
    except:
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


def extract_nanopore_for_sample(args, sample, filtered_reads, raw_nanopore_fasta):
    """
    Function aligns the Illumina reads to long reads and creates a new FASTA file of reads that were aligned to
    :param args:
    :param sample:
    :param filtered_reads:
    :param raw_nanopore_fasta:
    :return: Name of the output FASTA file
    """
    logging.info("Aligning Illumina reads to long reads... ")

    select_nanopore = sample.output_dir + os.sep + sample.id + "_select_nanopore.fasta"
    check_index(args.nanopore_reads, args.executables["bwa"])

    filtered_forward = ""
    filtered_reverse = ""
    for fastq in filtered_reads:
        if fastq.find("BackboneFiltered_R_1"):
            filtered_forward = fastq
        if fastq.find("BackboneFiltered_R_2"):
            filtered_reverse = fastq
        else:
            raise AssertionError("Unknown file found in filtered reads list: " + fastq)

    align_command = [args.executables["bwa"], "mem", "-t", args.threads, args.nanopore_reads, filtered_forward]
    if args.reverse:
        align_command.append(filtered_reverse)
    align_command.append("1>")
    sam_file = sample.output_dir + sample.id + "_nanopore_hits.sam"
    align_command.append(sam_file)
    align_command += ["2>", "/dev/null"]
    p_mem = subprocess.Popen(' '.join(align_command), shell=True, preexec_fn=os.setsid)
    p_mem.wait()
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


def correct_nanopore(args, sample):
    """
    A function to correct errors in the nanoopore reads using short, accurate Illumina reads and proovread
    :param sample: Miffed object with information of current sample
    :param args: command-line arguments list
    :return: Name of the corrected and trimmed nanopore reads
    """
    logging.info("Correcting errors in " + args.nanopore_reads + " using proovread... ")

    proovread_prefix = sample.output_dir + sample.id + "_proovread."
    proovread_command = [args.executables["proovread"], "--threads", args.threads]
    proovread_command += ["--long-reads=" + args.nanopore_reads]
    proovread_command += ["--short-reads=" + sample.forward_reads]
    if args.reverse:
        proovread_command += ["--short-reads=" + sample.reverse_reads]
    proovread_command += ["-p", sample.output_dir + os.sep + "proovread"]  # prefix to output files
    try:
        correct_stderr = open(proovread_prefix + "stderr", 'w')
    except:
        raise IOError("ERROR: cannot open file: " + proovread_prefix + "stderr")
    try:
        correct_stdout = open(proovread_prefix + "stdout", 'w')
    except:
        raise IOError("ERROR: cannot open file: " + proovread_prefix + "stdout")

    p_correct = subprocess.Popen(' '.join(proovread_command), shell=True, preexec_fn=os.setsid,
                                 stderr=correct_stderr, stdout=correct_stdout)
    p_correct.wait()
    correct_stderr.close()
    correct_stdout.close()
    logging.info("done.\n")
    return sample.output_dir + os.sep + "proovread" + os.sep + "proovread.trimmed.fa"


def quality_trimming(sample, args, filtered_reads):
    """
    Wrapper for trimmomatic
    :param sample: Miffed object with information of current sample
    :param args: command-line arguments list
    :param filtered_reads: list of background-filtered fastq files
    :return: list of quality-trimmed fastq files
    """
    logging.info("Trimming reads... ")
    trimmomatic_command = ["java", "-jar", args.executables["trimmomatic"], "PE", "-threads", args.threads]
    trimmomatic_command += filtered_reads
    trim_prefix = sample.output_dir + os.sep + sample.id + "_trim_"
    trimmomatic_outputs = [trim_prefix + "pe.1.fq", trim_prefix + "se.1.fq",
                           trim_prefix + "pe.2.fq", trim_prefix + "se.2.fq"]
    trimmomatic_command += trimmomatic_outputs
    # (/home/cmorganlang/bin/Trimmomatic-0.35/adapters/TruSeq3-PE.fa)
    adapters = "ILLUMINACLIP:" + args.adaptor_trim + "TruSeq3-PE.fa:2:3:10"
    trimmomatic_command.append(adapters)
    trimmomatic_command += ["LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"]
    try:
        trim_stderr = open(trim_prefix + "stderr.txt", 'w')
    except IOError:
        logging.error("Cannot open file: " + trim_prefix + "stderr.txt for writing.\n")
        sys.exit(3)
    try:
        trim_stdout = open(trim_prefix + "stdout.txt", 'w')
    except IOError:
        logging.error("Cannot open file: " + trim_prefix + "stdout.txt for writing.\n")
        sys.exit(3)
    p_trim = subprocess.Popen(' '.join(trimmomatic_command), shell=True, preexec_fn=os.setsid,
                              stderr=trim_stderr, stdout=trim_stdout)
    p_trim.wait()
    trim_stderr.close()
    trim_stdout.close()
    logging.info("done.\n")
    return trimmomatic_outputs


def determine_k_values(reads, assembler):
    k_min = 71
    if assembler == "spades":
        k_max = 127
    else:
        k_max = 241
    max_read_length = 0
    read_lengths = list()
    # Sample the first paired-end FASTQ file
    test_fastq = reads["pe_trimmed"][0]

    try:
        fq_handler = open(test_fastq)
    except IOError:
        logging.error("Unable to open FASTQ file '" + test_fastq + "' for reading.\n")
        sys.exit(3)

    x = 0
    for read in readfq(fq_handler):
        _, seq, _ = read
        read_length = len(seq)
        read_lengths.append(read_length)
        if read_length > max_read_length:
            max_read_length = read_length
        if x > 1E5:
            break
        x += 1
    fq_handler.close()

    if max_read_length < k_max:
        if max_read_length % 2 == 0:
            k_max = max_read_length - 3
        else:
            k_max = max_read_length - 2
        if max_read_length < k_min:
            k_min = k_max
    return k_min, k_max


def assemble_fosmids(sample: Sample, args, assembly_reads: dict, k_min: int, k_max: int, min_count: int):
    """
    Wrapper function for the assembly process - multi-sized de Bruijn graph based assembler for metagenomes
    :param sample: Miffed object with information of current sample
    :param args: command-line arguments list
    :param assembly_reads: Dictionary containing paired-end and single-end input reads for MEGAHIT
    :param k_min: 71; painstakingly determined to be optimal for fosmids sequenced with Illumina
    :param k_max: 241; painstakingly determined to be optimal for fosmids sequenced with Illumina
    :param min_count: The minimum k-mer abundance to be used by megahit for building the succinct DBG
    :return:
    """
    logging.info("Assembling sequences using " + args.assembler + "\n" +
                 "Parameters:\n--k-min = " + str(k_min) + "\t--k-max = " + str(k_max) +
                 "\t--k-step = 10\t--min-count = " + str(min_count) + "\n")

    if args.assembler == "spades":
        # The following is for assembling with SPAdes:
        forward = ""
        reverse = ""
        for fastq in assembly_reads["pe_trimmed"]:
            if re.search(r'pe.1.fq$', fastq):
                forward = fastq
            elif re.search(r'pe.2.fq$', fastq):
                reverse = fastq
            else:
                logging.error("Unrecognized FASTQ file name: " + fastq + "!\n")
                sys.exit(3)
        if not forward or not reverse:
            logging.error("Unable to find the paired-end FASTQ files for assembling" + sample.id + ".\n")
            sys.exit(3)

        spades_command = ["spades.py"]
        spades_command += ["-1", forward]
        spades_command += ["-2", reverse]
        spades_command += ["-s", assembly_reads["singletons_fq"]]
        spades_command += ["--careful"]
        spades_command += ["--memory", str(10)]
        spades_command += ["--threads", str(args.threads)]
        spades_command += ["--cov-cutoff", str(min_count)]
        spades_command += ["-k", ','.join([str(k) for k in range(k_min, k_max, 10)])]
        spades_command += ["-o", sample.output_dir + "assembly"]
        p_spades = subprocess.Popen(' '.join(spades_command), shell=True, preexec_fn=os.setsid)
        p_spades.wait()
        if p_spades.returncode != 0:
            logging.error("SPAdes did not complete successfully! Continuing on to the next sample.\n")

        contigs_file = sample.output_dir + "assembly" + os.sep + "contigs.fasta"
        asm_log = sample.output_dir + "assembly" + os.sep + "spades.log"
        scaffolds = sample.output_dir + "assembly" + os.sep + "scaffolds.fasta"

    elif args.assembler == "megahit":
        if not os.path.isfile(assembly_reads["se_fa"]):
            logging.error("FASTA file containing trimmed orphaned reads could not be located!\n")
            sys.exit(3)
        if not os.path.isfile(assembly_reads["pe_fa"]):
            logging.error("FASTA file containing trimmed paired-end reads could not be located!\n")
            sys.exit(3)

        megahit_command = [args.executables["megahit"]]
        megahit_command += ["--12", assembly_reads["pe_fa"]]
        megahit_command += ["--read", assembly_reads["se_fa"]]
        megahit_command += ["--k-min", str(k_min), "--k-max", str(k_max)]
        megahit_command += ["--min-count", str(min_count)]
        megahit_command += ["--k-step", str(10)]
        megahit_command += ["--memory", str(0.25)]
        megahit_output_dir = sample.output_dir + "assembly"
        megahit_command += ["--out-dir", megahit_output_dir]
        megahit_command += ["--num-cpu-threads", args.threads]
        megahit_command.append("--verbose")
        megahit_command.append("--no-mercy")  # Recommended for high-coverage datasets
        megahit_command += ["1>", "/dev/null", "2>", "/dev/null"]
        # Run the megahit command using subprocess
        p_megahit = subprocess.Popen(' '.join(megahit_command), shell=True, preexec_fn=os.setsid)
        p_megahit.wait()
        if p_megahit.returncode != 0:
            logging.error("MEGAHIT did not complete successfully!\n")

        contigs_file = sample.output_dir + "assembly" + os.sep + "final.contigs.fa"
        asm_log = sample.output_dir + "assembly" + os.sep + "log"
        scaffolds = ""
    else:
        logging.error("Unknown assembly software '" + args.assembler + "' requested.\n")
        sys.exit(3)

    logging.info("Cleaning up assembly outputs... ")
    if not os.path.isfile(contigs_file):
        logging.error(args.assembler + " assembly was not created! Check " +
                      sample.output_dir + "assembly" + os.sep + "log for an error.")
        sys.exit(3)
    shutil.move(contigs_file,
                sample.output_dir + os.sep + sample.id + "_contigs.fasta")
    shutil.move(asm_log, sample.output_dir + os.sep + args.assembler + "_log.txt")
    shutil.rmtree(sample.output_dir + "assembly" + os.sep)
    if scaffolds:
        # Only available for SPAdes
        shutil.move(scaffolds,
                    sample.output_dir + os.sep + sample.id + "_scaffolds.fasta")
    logging.info("done.\n")

    return


def prep_reads_for_assembly(sample, args, trimmed_reads):
    logging.info("Preparing quality FASTQ files for assembly... ")
    reads = dict()
    reads["pe_trimmed"] = []
    reads["se_trimmed"] = []
    reads["singletons_fq"] = sample.output_dir + "singletons.fastq"

    for fastq in trimmed_reads:
        if re.search(r'pe.[1-2].fq$', fastq):
            reads["pe_trimmed"].append(fastq)
        else:
            reads["se_trimmed"].append(fastq)

    try:
        orphans_handler = open(reads["singletons_fq"], 'w')
    except IOError:
        logging.error("Unable to open " + reads["singletons_fq"] + " for writing.\n")
        sys.exit(3)
    for unpaired_file in reads["se_trimmed"]:
        with open(unpaired_file) as orphan_fq:
            for line in orphan_fq:
                orphans_handler.write(line)
    orphans_handler.close()

    if args.assembler == "megahit":
        # MEGAHIT requires FASTA files so interleave and convert to FASTA using fq2fa
        reads["pe_fa"] = sample.output_dir + sample.id + "_paired.fa"
        reads["se_fa"] = sample.output_dir + sample.id + "_unpaired.fa"
        paired_merge = [args.executables["fq2fa"], "--merge"]
        paired_merge += reads["pe_trimmed"]
        paired_merge.append(reads["pe_fa"])
        p_pemerge = subprocess.Popen(' '.join(paired_merge), shell=True, preexec_fn=os.setsid)
        p_pemerge.wait()

        unpaired_fq2fa = [args.executables["fq2fa"], reads["singletons_fq"], reads["se_fa"]]
        p_sefq2fa = subprocess.Popen(' '.join(unpaired_fq2fa), shell=True, preexec_fn=os.setsid)
        p_sefq2fa.wait()

    logging.info("done.\n")
    return reads


def deinterleave_fastq(fastq_file, output_dir=""):
    try:
        fq_in_handler = open(fastq_file, 'r')
    except IOError:
        logging.error("Unable to open " + fastq_file + " for reading.\n")
        sys.exit(3)

    logging.info("De-interleaving forward and reverse reads in " + fastq_file + "... ")
    if not output_dir:
        output_dir = os.path.dirname(fastq_file)
    prefix = '.'.join(os.path.basename(fastq_file).split('.')[:-1])
    fwd_fq = output_dir + os.sep + prefix + "_R1.fastq"
    rev_fq = output_dir + os.sep + prefix + "_R2.fastq"

    try:
        f_handler = open(fwd_fq, 'w')
        r_handler = open(rev_fq, 'w')
    except IOError:
        logging.error("Unable to open deinterleaved FASTQ files for writing in " + output_dir + "\n")
        sys.exit(3)

    f_string = ""
    r_string = ""
    acc = 0
    for name, seq, qual in readfq(fq_in_handler):
        if acc % 2:
            r_string += "\n".join(["@" + name, seq, '+', qual]) + "\n"
        else:
            f_string += "\n".join(["@" + name, seq, '+', qual]) + "\n"

        acc += 1
        if acc % 1E6 == 0:
            f_handler.write(f_string)
            r_handler.write(r_string)

    # Close up shop
    f_handler.write(f_string)
    f_handler.close()
    r_handler.write(r_string)
    r_handler.close()
    logging.info("done.\n")

    return fwd_fq, rev_fq


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
        logging.error("Unable to locate fastq files. Must end in either 'fastq' or 'fq'\n")
        sys.exit(3)
    regex_sample = re.compile(sample_id)
    for fastq in forward_reads:
        if regex_sample.search(os.path.basename(fastq)) and raw_reads["forward"] == "":
            raw_reads["forward"] = os.path.join(os.getcwd(), fastq)
        # Ensure there is a single forward fastq file for sample_id
        elif regex_sample.search(os.path.basename(fastq)) and raw_reads["forward"] != "":
            logging.error("More than 2 fastq files match sample ID string in reads directory! " +
                          "Please concatenate those files from the same library and try again.\n")
            sys.exit(3)
        else:
            # FASTQ file is not of the current sample
            pass
    if args.reverse:
        reverse_reads = glob.glob(args.reverse + os.sep + "*fastq")
        reverse_reads += glob.glob(args.reverse + os.sep + "*fq")
        for fastq in reverse_reads:
            if re.search(sample_id, fastq) and raw_reads["reverse"] == "":
                fastq = os.path.join(os.getcwd(), fastq)
                raw_reads["reverse"] = fastq
            # Ensure there is a single reverse fastq file for sample_id
            elif regex_sample.search(os.path.basename(fastq)) and raw_reads["reverse"] != "":
                logging.error("More than 2 fastq files match sample ID string in reverse directory! " +
                              "Please concatenate those files from the same library and try again.\n")
                sys.exit(3)

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


def get_assemble_input(sample):
    response = input("Assembled fosmids for " + sample.id + " already exists. Do you want to re-assemble? [y|n]")
    if response == "y":
        assemble = True
        logging.info("Re-assembling " + str(sample.id) + "\n")
    elif response == "n":
        assemble = False
        logging.info("Using " + str(sample.assembled_fosmids) + "\n")
    else:
        logging.info("Unknown response" + "\n")
        assemble = get_overwrite_input(sample.id)
    return assemble


def parse_miffed(args):
    """
    Finds the primary project name, submitter, sequencing platform and other information to help with processing
    and storage of the inputs/outputs.
    :param args: parsed command-line arguments from get_options()
    :return:
    """
    miffed = open(args.miffed, 'r')
    line = miffed.readline()
    libraries = list()
    header_re = re.compile(r"Sample.*Project.*Human selector.*Number of fosmids")

    while line:
        if header_re.search(line):
            pass
        elif not line[0] == '#':
            if ',' not in line:
                logging.error("MIFFED input file doesn't seem to be in csv format!\n")
                sys.exit(3)
            fields = line.strip().split(',')
            if len(fields) < 13:
                raise AssertionError(str(len(fields)) + " fields in " + args.miffed +
                                     " when at least 13 are expected\n" + line)
            miffed_entry = Miffed(fields[0])
            miffed_entry.populate_info(fields)
            miffed_entry.output_dir = os.sep.join([args.fabfos_path, miffed_entry.project, miffed_entry.id]) + os.sep
            miffed_entry.assembled_fosmids = miffed_entry.output_dir + os.sep + miffed_entry.id + "_contigs.fasta"
            complete = miffed_entry.ensure_completeness(args.miffed)
            if complete is False:
                miffed_entry.exclude = True
            if args.nanopore_reads:
                miffed_entry.nanopore = True
            if args.skip_correction is True:
                miffed_entry.error_correction = False

            # Check to ensure the output directory exist or ensure it is empty:
            if os.path.isdir(miffed_entry.output_dir):
                files = glob.glob(miffed_entry.output_dir + "*")
                if len(files) != 0:
                    miffed_entry.overwrite = get_overwrite_input(miffed_entry.id)
                    if miffed_entry.overwrite:
                        try:
                            shutil.rmtree(miffed_entry.output_dir)
                            os.makedirs(miffed_entry.output_dir)
                        except Exception as e:
                            logging.error(str(e) + "\n")
                            raise
                    else:
                        # Check to see if sample should be excluded entirely or not
                        miffed_entry.exclude = get_exclude_input(miffed_entry.id)

                        if not miffed_entry.exclude:
                            if os.path.isfile(miffed_entry.assembled_fosmids):
                                miffed_entry.assemble = get_assemble_input(miffed_entry)
            else:
                os.makedirs(miffed_entry.output_dir)

            project_metadata_file = args.fabfos_path + os.sep + miffed_entry.project + os.sep +\
                                    "FabFos_" + miffed_entry.project + "_metadata.tsv"
            if not os.path.isfile(project_metadata_file):
                project_metadata = open(project_metadata_file, 'w')
                project_metadata.write(args.metadata_header)
                project_metadata.close()
            libraries.append(miffed_entry)
        line = miffed.readline()
    miffed.close()
    return libraries, args


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
    """
    Function to count the number of reads in all FASTQ files in file_list
    :param file_list: list of FASTQ files
    :return: integer representing the number of reads in all FASTQ files provided
    """
    num_reads = 0
    for fastq in file_list:
        if not fastq:
            continue
        fq_in = open(fastq, 'r')
        for read_dat in readfq(fq_in):
            num_reads += 1
        fq_in.close()
    return num_reads


def map_ends(args, sample):
    """
    Function for using blastn to align fosmid ends in a file to find the well of each contig
    """
    # Index the contigs
    sample_prefix = sample.output_dir + os.sep + sample.id
    blastdb_command = [args.executables["makeblastdb"]]
    blastdb_command += ["-in", sample_prefix + "_contigs.fasta"]
    blastdb_command += ["-dbtype", "nucl"]
    blastdb_command += ["-out", sample_prefix]
    blastdb_command += ["1>/dev/null", "2>/dev/null"]
    p_makeblastdb = subprocess.Popen(" ".join(blastdb_command), shell=True, preexec_fn=os.setsid)
    p_makeblastdb.wait()

    logging.info("Aligning fosmid ends to assembly... ", "out")
    blastn_command = [args.executables["blastn"]]
    blastn_command += ["-db", sample_prefix]
    blastn_command += ["-outfmt", "\"6",
                       "sseqid", "slen", "qseqid", "pident", "length", "sstrand",
                       "sstart", "send", "bitscore", "qstart", "qend", "\""]
    blastn_command += ["-query", args.ends]
    blastn_command += ["-perc_identity", str(95)]
    blastn_command += ["-out", sample_prefix + "_endsMapped.tbl"]
    blastn_command += ["1>/dev/null", "2>/dev/null"]
    p_align_ends = subprocess.Popen(" ".join(blastn_command), shell=True, preexec_fn=os.setsid)
    p_align_ends.wait()
    logging.info("done\n")

    return


def get_fosmid_end_name(string):
    """
    This function is very unstable and currently needs to be manually changed for every FabFos run
    """
    # PPSLIBM-08-H17_F.ab1
    full_name = string.split('.')[0]
    prefix = full_name.split('-')[0]
    # print "string= ", string, "\tfull= ", full_name, "\tprefix= ", prefix
    name, direction = full_name.split('_')
    # name = prefix + '.' + name
    return name


def parse_end_alignments(sample, fosmid_assembly, ends_stats):
    """
    Function parses BLAST alignments of fosmid ends mapped to the MEGAHIT assemblies
    :param sample: Miffed object with information of current sample
    :param fosmid_assembly: dictionary with headers as keys and sequences for values
    :param ends_stats: dictionary containing information on the fosmid ends
    :return: dictionary of contig names (keys) and list of EndAlignment objects (values) or empty list
    """
    ends_mapping = dict()
    ends_stats["aligned"] = set()
    ends_stats["unaligned"] = set()
    blast_output = sample.output_dir + os.sep + sample.id + "_endsMapped.tbl"
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
        logging.warning("Forward strand did not align in single-pair clone " + clone + "\n")
    if not reverse_aligned:
        logging.warning("Reverse strand did not align in single-pair clone " + clone + "\n")
    if len(positions_dict["minus_strand"]) == 0 or len(positions_dict["plus_strand"]) == 0:
        logging.warning("Forward and reverse mate-pairs of " + clone + " aligned to the same strand of " + contig +
                        "!\nThis sequence should not be trusted due to assembly or library preparation errors\n")

    positions = positions_dict["plus_strand"] + positions_dict["minus_strand"]
    start = min(positions)
    end = max(positions)
    return start, end


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


def prune_and_scaffold_fosmids(sample, clone_map, multi_fosmid_map):
    fragments_tsv = sample.output_dir + os.sep + "potential_fosmid_fragments.tsv"

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


def write_fosmid_assignments(sample, clone_map):
    coord_name = sample.output_dir + os.sep + sample.id + "_fosmid_coordinates.tsv"
    fasta_name = sample.output_dir + os.sep + sample.id + "_formatted.fasta"
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
    try:
        fasta_seqs = open(fasta, 'r')
    except IOError:
        logging.error("Cannot open FASTA file (" + fasta + ") for reading!")
        sys.exit(3)

    line = fasta_seqs.readline()
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
        line = fasta_seqs.readline()
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


def get_assembly_stats(sample, args, fosmid_assembly):
    assembly_stats = dict()
    assembly_stats["Assembler"] = get_assembler_version(args.assembler)
    sample_prefix = sample.output_dir + os.sep + sample.id
    nx_command = [args.executables["getNx"]]
    nx_command += ["-i", sample_prefix + "_contigs.fasta"]
    nx_command += ["-o", sample_prefix + "_nx.csv"]
    nx_command += ["-m", str(2000)]
    nx_command.append(">" + sample.output_dir + os.sep + sample.id + "_nx.txt")

    p_nx = subprocess.Popen(' '.join(nx_command), shell=True, preexec_fn=os.setsid)
    p_nx.wait()
    if p_nx.returncode != 0:
        logging.error("getNx did not finish successfully!\n")
        sys.exit(3)

    try:
        nx_stats = open(sample_prefix + "_nx.csv", 'r')
    except IOError:
        logging.error("getNx did not finish successfully: unable to open " + sample_prefix + "_nx.csv for reading!\n")
        sys.exit(3)

    nx_stats.readline()
    line = nx_stats.readline()
    while line:
        line = line.strip()
        proportion, n = line.split(',')
        if float(proportion) == 0.5:
            logging.info("N50 = " + str(n) + "\n")
            assembly_stats["N50"] = str(n)
        if float(proportion) == 0.9:
            assembly_stats["N90"] = str(n)
        if float(proportion) == 0:
            logging.info("Longest contig = " + str(n) + "bp\t")
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

    os.remove(sample.output_dir + os.sep + sample.id + "_nx.txt")

    return assembly_stats


def update_metadata(metadata_file, sample_id, library_metadata, read_stats, assembly_stats, ends_stats):
    try:
        metadata = open(metadata_file, 'a')
    except IOError:
        logging.error("Unable to open " + metadata_file + " to append library metadata!\n")
        sys.exit(3)

    metadata.write(sample_id + "\t")
    metadata.write(library_metadata.project + "\t")
    metadata.write(library_metadata.selector + "\t")
    metadata.write(library_metadata.vector + "\t")
    metadata.write(library_metadata.screen + "\t")
    metadata.write(library_metadata.selection + "\t")
    metadata.write(str(library_metadata.num_fosmids) + "\t")
    metadata.write(library_metadata.seq_submission_date + "\t")
    metadata.write(library_metadata.glycerol_plate + "\t")
    metadata.write(library_metadata.seq_center + "\t")
    metadata.write(library_metadata.seq_type + "\t")
    metadata.write(str(library_metadata.read_length) + "\t")
    metadata.write(library_metadata.instrument + "\t")
    metadata.write(strftime("%Y-%m-%d") + "\t")
    metadata.write(read_stats["num_raw_reads"] + "\t")
    metadata.write(read_stats["percent_off_target"] + "\t")
    metadata.write(str(read_stats["number_trimmed_reads"]) + "\t")
    metadata.write(assembly_stats["Assembler"] + "\t")
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


def write_fosmid_end_failures(sample, ends_stats):
    """
    Writes the names of fosmid ends that could not be aligned and failed
    :param sample: Miffed object with information of current sample
    :param ends_stats:
    :return:
    """
    sample_prefix = sample.output_dir + os.sep + sample.id
    fosmid_end_failures = sample_prefix + "_fosmid_end_failures.tsv"

    try:
        failure_file = open(fosmid_end_failures, 'w')
    except IOError:
        logging.error("Cannot open file: " + fosmid_end_failures)
        sys.exit(3)
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


def write_unique_fosmid_ends_to_bulk(args):
    logging.info("Writing new fosmid end sequences to the legacy fosmid end file... ")
    new_fosmid_ends = read_fasta(args.ends)
    bulk_ends_file = args.fabfos_path + "FabFos_legacy_ends.fasta"
    bulk_ends = read_fasta(bulk_ends_file)

    try:
        legacy_ends = open(bulk_ends_file, 'a')
    except:
        raise IOError("Cannot open " + bulk_ends_file + " for appending!")

    for new_end in new_fosmid_ends:
        if new_end not in bulk_ends.keys():
            if len(new_fosmid_ends[new_end]) > 100:
                legacy_ends.write(new_end + "\n")
                legacy_ends.write(new_fosmid_ends[new_end] + "\n")

    legacy_ends.close()
    logging.info("done.\n")

    return


def align_nanopore_to_background(args, sample):
    logging.info("Aligning " + sample.nanopore + " to background sequences... ")

    # Make the BLAST database
    try:
        background_blastdb_stderr = open(sample.output_dir + "background_blastdb.stderr", 'w')
        background_blastdb_stdout = open(sample.output_dir + "background_blastdb.stdout", 'w')
    except IOError:
        logging.error("Unable to open " + sample.output_dir + "background_blastdb.stderr for writing")
        sys.exit(3)

    background_db = args.background + "_BLAST"
    blast_db_command = [args.executables["makeblastdb"]]
    blast_db_command += ["-dbtype", "nucl"]
    blast_db_command += ["-in", args.background]
    blast_db_command += ["-out", background_db]
    p_makeblastdb = subprocess.Popen(' '.join(blast_db_command), shell=True, preexec_fn=os.setsid,
                                     stderr=background_blastdb_stderr, stdout=background_blastdb_stdout)
    p_makeblastdb.wait()
    background_blastdb_stderr.close()
    background_blastdb_stdout.close()

    # Align the corrected reads to the BLAST database
    nanopore_background_alignments = sample.output_dir + os.sep + "nanopore_background_BLAST.tbl"
    blastn_command = [args.executables["blastn"]]
    blastn_command += ["-query", sample.nanopore]
    blastn_command += ["-db", background_db]
    blastn_command += ["-outfmt", "\"6",
                       "sseqid", "slen", "qseqid", "pident", "length", "sstrand",
                       "sstart", "send", "bitscore", "qstart", "qend", "\""]
    blastn_command += ["-out", nanopore_background_alignments]
    p_blastn = subprocess.Popen(' '.join(blastn_command), shell=True, preexec_fn=os.setsid)
    p_blastn.wait()

    # Remove the BLAST database
    blastdb_files = glob.glob('.'.join(args.background.split('.')[0:-1]) + "_BLAST.*")
    for db_file in blastdb_files:
        os.remove(db_file)

    logging.info("done.\n")
    return nanopore_background_alignments


def trim_background(sample, nanopore_background_alignments, read_stats):
    """
    Function
    :param nanopore_background_alignments:
    :param sample: Miffed object with information of current sample
    :param read_stats: A dictionary with values representing basic stats on number of reads and trimming
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
    except:
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
    read_stats["percent_off_target"] = str(truncated/int(read_stats["num_raw_reads"]))
    read_stats["number_trimmed_reads"] = str(truncated)

    return nanopore_reads, read_stats


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
    trimmed_nanopore_fasta = miffed_entry.output_dir + os.sep + miffed_entry.id + "_nanopore_trimmed.fasta"
    try:
        new_fasta = open(trimmed_nanopore_fasta, 'w')
    except:
        raise IOError("Unable to open " + trimmed_nanopore_fasta + " for writing!")

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


def determine_min_count(num_reads, num_fosmids, k_max):
    """
    Function to determine the best value for the minimum coverage of a k-mer to be included in assembly
    :param num_reads:
    :param num_fosmids:
    :param k_max:
    :return:
    """
    approx_coverage = (num_reads * k_max) / (int(num_fosmids) * 40000)
    sys.stdout.write("Approximate fosmid coverage = " + str(approx_coverage) + "\n")
    sys.stdout.flush()
    min_count = 10
    dist_tail = approx_coverage / 100
    if dist_tail > min_count:
        min_count = int(dist_tail)

    return min_count


def assemble_nanopore_reads(sample, args):
    """
    Wrapper function for launching canu assembler with the corrected nanopore reads
    :param sample:
    :param args:
    :return:
    """
    logging.info("Assembling error-corrected nanopore reads with canu... ", "out")

    canu_output = sample.output_dir + os.sep + "canu" + os.sep
    if os.path.isdir(canu_output):
        shutil.rmtree(canu_output)
    try:
        os.makedirs(canu_output)
    except:
        raise IOError("Unable to make " + canu_output)

    corrected_nanopore_fasta = sample.output_dir + os.sep + sample.id + "_nanopore_trimmed.fasta"
    canu_stderr = open(canu_output + "canu.stderr", 'w')
    canu_stdout = open(canu_output + "canu.stdout", 'w')

    genome_size = int(sample.num_fosmids) * 40

    canu_command = [args.executables["canu"]]
    canu_command += ["genomeSize=" + str(genome_size) + "k", "minThreads=" + str(args.threads)]
    if sample.error_correction is True:
        # TODO: Compare assemblies where proovread-corrected nanopore reads are just assembled or ran through all stages
        canu_command.append("-assemble")
    canu_command += ["-p", sample.id]
    canu_command += ["-d", canu_output]
    if args.skip_correction:
        canu_command += ["-nanopore-raw", corrected_nanopore_fasta]
    else:
        canu_command += ["-nanopore-corrected", corrected_nanopore_fasta]

    p_canu = subprocess.Popen(' '.join(canu_command), shell=True, preexec_fn=os.setsid,
                              stderr=canu_stderr, stdout=canu_stdout)
    p_canu.wait()

    shutil.move(canu_output + sample.id + ".contigs.fasta", sample.assembled_fosmids)

    logging.info("done.\n")

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
        ends_stats["Num_total"] = 0
        ends_stats["Num_failed"] = 0
        ends_stats["Num_aligned"] = 0
        ends_stats["Num_unaligned"] = 0

    # Sequence processing and filtering:
    for sample in libraries:
        read_stats = dict()
        read_stats["num_raw_reads"] = "NA"
        read_stats["percent_off_target"] = "NA"
        read_stats["number_trimmed_reads"] = "NA"
        if not sample.exclude:
            if sample.assemble:
                logging.info("Processing raw data for " + sample.id + "\n" +
                             "Outputs for " + sample.id + " will be found in " + sample.output_dir + "\n")
                raw_reads = find_raw_reads(args, sample.id)
                if args.interleaved:
                    raw_reads["forward"], raw_reads["reverse"] = deinterleave_fastq(raw_reads["forward"],
                                                                                    sample.output_dir)
                sample.forward_reads = raw_reads["forward"]
                sample.reverse_reads = raw_reads["reverse"]
                num_raw_reads = find_num_reads(raw_reads.values())
                if num_raw_reads == 0:
                    logging.error("No reads found for " + sample.id + "\n")
                    sys.exit(3)
                read_stats["num_raw_reads"] = str(num_raw_reads)
                logging.info("Number of raw reads = " + str(num_raw_reads) + "\n")
                filtered_reads, bam_file = filter_backbone(sample, args, raw_reads)
                clustering(bam_file, sample.output_dir, args.bb_length, args.cluster_length)

                num_filtered_reads = find_num_reads(filtered_reads)
                logging.info(str(num_raw_reads - num_filtered_reads) + " reads removed from background filtering (" +
                             str(((num_raw_reads - num_filtered_reads) * 100) / num_raw_reads) + "%).\n")
                read_stats["percent_off_target"] = str(
                    ((num_raw_reads - num_filtered_reads) * 100) /
                    num_raw_reads)
                if num_filtered_reads < 1600:
                    # The number of reads remaining is too low for assembly (< 20X for a single fosmid)
                    logging.warning("Number of reads remaining will provide less than 20X coverage for a single fosmid"
                                    " - skipping this sample\n")
                    continue
                trimmed_reads = quality_trimming(sample, args, filtered_reads)
                num_reads_assembly = find_num_reads(trimmed_reads)
                read_stats["number_trimmed_reads"] = str(
                    (num_filtered_reads - num_reads_assembly) * 100 /
                    num_reads_assembly)
                assembly_reads = prep_reads_for_assembly(sample, args, trimmed_reads)

                if sample.nanopore:
                    raw_nanopore_fasta = read_fasta(args.nanopore_reads)
                    read_stats["num_raw_reads"] = str(len(raw_nanopore_fasta.keys()))
                    if sample.error_correction is False:
                        # Skip the error correction and filter the nanopore reads that don't align to Illumina reads
                        # Set nanopore_reads to the filtered raw reads
                        sample.nanopore = extract_nanopore_for_sample(args, sample, filtered_reads, raw_nanopore_fasta)
                    else:
                        sample.nanopore = correct_nanopore(args, sample)
                    # Align the corrected reads to the trim_sequences.fasta file using LAST
                    nanopore_background_alignments = align_nanopore_to_background(args, sample)
                    # Trim the background sequences and reads shorter than 1000bp after removing background
                    nanopore_reads, read_stats = trim_background(sample, nanopore_background_alignments, read_stats)
                    write_trimmed_reads(sample, nanopore_reads)
                    # Assemble nanopore reads
                    assemble_nanopore_reads(sample, args)

                else:
                    # TODO: Include an optional minimus2 module to further assemble the contigs
                    k_min, k_max, = determine_k_values(assembly_reads, args.assembler)
                    min_count = determine_min_count(num_reads_assembly, sample.num_fosmids, k_max)
                    assemble_fosmids(sample, args, assembly_reads, k_min, k_max, min_count)
                    clean_intermediates(sample)

            fosmid_assembly = read_fasta(sample.assembled_fosmids)
            assembly_stats = get_assembly_stats(sample, args, fosmid_assembly)
            # For mapping fosmid ends:
            if args.ends:
                map_ends(args, sample)
                ends_mapping, ends_stats = parse_end_alignments(sample, fosmid_assembly, ends_stats)
                clone_map, ends_stats, multi_fosmid_map = assign_clones(ends_mapping, ends_stats, fosmid_assembly)
                clone_map = prune_and_scaffold_fosmids(sample, clone_map, multi_fosmid_map)
                write_fosmid_assignments(sample, clone_map)
                write_fosmid_end_failures(sample, ends_stats)
                write_unique_fosmid_ends_to_bulk(args)
                # TODO: Remove blast database for ends
                # TODO: include a visualization to show where the fosmid ends mapped on each contig.
            project_metadata_file = args.fabfos_path + os.sep + sample.project + os.sep + \
                                    "FabFos_" + sample.project + "_metadata.tsv"
            update_metadata(project_metadata_file, sample.id, sample, read_stats, assembly_stats, ends_stats)
            update_metadata(args.master_metadata, sample.id, sample, read_stats, assembly_stats, ends_stats)


main()
