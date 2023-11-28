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

import json
import os, sys
from pathlib import Path
import argparse
import inspect
from dataclasses import dataclass
import multiprocessing
import importlib

from .models import Assembly, BackgroundGenome, EndSequences, ReadsManifest, VectorBackbone
from .utils import NAME, USER, VERSION, ENTRY_POINTS, MODULE_ROOT

CLI_ENTRY = ENTRY_POINTS[0]
    
class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help(sys.stderr)
        self.exit(2, '\n%s: error: %s\n' % (self.prog, message))


class CommandLineInterface:
    def _get_fn_name(self):
        return inspect.stack()[1][3]

    def run(self, raw_args):
        parser = ArgumentParser(
            prog = f'{CLI_ENTRY} {self._get_fn_name()}',
        )

        # the arguments here are tightly coupled to models.py, sorry
        DEFAULT_ASSEMBLY_MODES = Assembly.CHOICES[:2]
        paths = parser.add_argument_group(title="main")
        paths.add_argument("-1", "--forward", metavar="FASTQ", nargs='*', required=False, default=[],
            help="forward paired-end reads")
        paths.add_argument("-2", "--reverse", metavar="FASTQ", nargs='*', required=False, default=[],
            help="reverse paired-end reads in the same order")
        paths.add_argument("-s", "--single", metavar="FASTQ", nargs='*', required=False, default=[],
            help="single-end reads")
        paths.add_argument("-i", "--interleaved", metavar="FASTQ", nargs='*', required=False, default=[],
            help="interleaved reads")
        paths.add_argument("-o", "--output", metavar="PATH", required=True,
            help="path to output folder, will be created if non-existent")
        
        fos = parser.add_argument_group(title="fosmid pool specific")
        fos.add_argument("--no_trim", action="store_true", default=False, required=False,
            help="skip read trimming with trimmomatic")
        fos.add_argument("-b", "--background", metavar="FASTA", required=False,
            help="host background to filter out")
        fos.add_argument("--endf", metavar="FASTA", nargs='*', required=False,
            help="sanger end sequences")
        fos.add_argument("--endr", metavar="FASTA", nargs='*', required=False,
            help="sanger end sequences from the other end, IDs must match those from --endf")
        fos.add_argument("--end_regex", metavar="STR", required=False,
            help="regex for getting ID of end seq., default: \"%s\", ex. \"\\w+_\\d+\" would get ABC_123 from ABC_123_FW" % EndSequences.DEFAULT_REGEX)
        fos.add_argument("--vector", metavar="FASTA", required=False,
            help="the vector backbone sequence for pool size estimation")

        # "options" group
        parser.add_argument("-a", "--assemblies", nargs='*', required=False, default=[],
            help=f"pre-assembled contigs or assembly modes to use, pick any combination of {Assembly.CHOICES}, default:{DEFAULT_ASSEMBLY_MODES}")
        parser.add_argument("--min_length", metavar="INT", required=False, default=1000,
            help="min contig length to accept, default=1000")
        parser.add_argument("--overwrite", action="store_true", default=False, required=False,
            help="overwrite previous output, if given same output path")
        parser.add_argument("-t", "--threads", metavar="INT", type=int,
            help="threads, default:ALL", default=multiprocessing.cpu_count())
        parser.add_argument("--mock", action="store_true", default=False, required=False,
            help="dry run snakemake")
        parser.add_argument("--snakemake", nargs='*', required=False, default=[],
            help="additional snakemake cli args in the form of KEY=VALUE or KEY (no leading dashes)")
        args = parser.parse_args(raw_args)

        #########################
        # verify & parse inputs
        #########################
        input_error = False
        _printed = False
        def _error(message: str):
            nonlocal input_error, _printed
            if not _printed:
                parser.print_help()
                print()
                _printed = True
            print(f"Invalid input: {message}")
            input_error = True

        output = Path(args.output).absolute()
        logs = output.joinpath("logs")
        if not output.exists(): os.makedirs(output)
        if not logs.exists(): os.makedirs(logs)
        with open(output.joinpath("params.json"), "w") as j:
            d = args.__dict__|dict(
                current_directory=os.getcwd(),
            )
            for k in list(d):
                if isinstance(d[k], list) and len(d[k]) == 0: del d[k]
                elif d[k] is None: del d[k]
            json.dump(d, j, indent=4)

        input_models = {}
        for model_class in [
            ReadsManifest, BackgroundGenome, Assembly, EndSequences, VectorBackbone,
        ]:
            input_models[model_class] = model_class.Parse(args, output, _error)
        has_reads = len([r for g in input_models[ReadsManifest].AllReads() for r in g])>0
        selected_modes = len(input_models[Assembly].modes)>0
        given_assemblies = len(input_models[Assembly].given)>0
        if not has_reads and selected_modes:
            _error(f"selected assembly modes without giving reads")
        if not has_reads and not given_assemblies:
            _error(f"must provide reads, previously assembled contigs, or both")
        if has_reads and not selected_modes:
            input_models[Assembly].modes = DEFAULT_ASSEMBLY_MODES
            input_models[Assembly].Save(output.joinpath(Assembly.ARG_FILE))

        smk_args = ["--latency-wait 0"]
        for a in args.snakemake:
            if "=" in a:
                toks = a.split("=")
                pa = f"--{toks[0]} {'='.join(toks[1:])}"
            else:
                pa = f"--{a}"
            smk_args.append(pa)

        #########################
        # run snakemake
        #########################
        if input_error: return
        smk_log = logs.joinpath("snakemake")
        link_log = "" if smk_log.exists() else f'ln -s ../.snakemake/log {logs.joinpath("snakemake")}'
        params = dict(
            src=MODULE_ROOT,
            log=logs,
            threads=args.threads,
        )

        params_str = ' '.join(f"{k}={v}" for k, v in params.items())
        cache = output.joinpath("internals/temp_cache")
        cmd = f"""\
            {link_log}
            mkdir -p {cache}
            export XDG_CACHE_HOME={cache}
            snakemake -s {MODULE_ROOT.joinpath('main.smk')} -d {output} --rerun-incomplete \
                {'--forceall' if args.overwrite else ''} \
                {' '.join(smk_args)} \
                {"-n" if args.mock else ""} \
                --config {params_str} \
                --keep-going --keep-incomplete --cores {args.threads}
        """
        os.system(cmd)

    def api(self, raw_args=None):
        parser = ArgumentParser(
            prog = f'{CLI_ENTRY} {self._get_fn_name()}',
            description=f"Snakemake uses this to call the python script for each step"
        )

        parser.add_argument("--step", required=True)
        parser.add_argument("--args", nargs='*', required=False, default=[])
        args = parser.parse_args(raw_args)

        mo = importlib.import_module(name=f".steps.{args.step}", package=NAME)
        try:
            mo.Procedure(args.args)
        except KeyboardInterrupt:
            exit()

    def help(self, args=None):
        help = [
            f"{NAME} v{VERSION}",
            f"https://github.com/{USER}/{NAME}",
            f"",
            f"Syntax: {CLI_ENTRY} COMMAND [OPTIONS]",
            f"",
            f"Where COMMAND is one of:",
        ]+[f"- {k}" for k in COMMANDS]+[
            f"",
            f"for additional help, use:",
            f"{CLI_ENTRY} COMMAND -h/--help",
        ]
        help = "\n".join(help)
        print(help)
COMMANDS = {k:v for k, v in CommandLineInterface.__dict__.items() if k[0]!="_"}

def main():
    cli = CommandLineInterface()
    if len(sys.argv) <= 1:
        cli.help()
        return

    COMMANDS.get(# calls command function with args
        sys.argv[1], 
        CommandLineInterface.help # default
    )(cli, sys.argv[2:]) # cli is instance of "self"

if __name__ == "__main__":
    main()
