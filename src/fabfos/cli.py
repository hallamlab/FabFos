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


import os, sys
from pathlib import Path
import argparse
import inspect
from dataclasses import dataclass
import multiprocessing
import shutil
import importlib

from .models import ReadsManifest
from .constants import ORIGINAL_READS, SKIP_FILTER
from .utils import NAME, VERSION, ENTRY_POINTS, MODULE_ROOT

CLI_ENTRY = ENTRY_POINTS[0]

def _line():
    try:
        width = os.get_terminal_size().columns
    except:
        width = 32
    return "="*width
    
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
        ASSEMBLERS = "megahit, spades_meta, spades_isolate, spades_sc,".split(", ")

        paths = parser.add_argument_group(title="paths")
        paths.add_argument("-1", "--forward", metavar="FASTQ", nargs='*', required=False, default=[],
            help="forward paired-end reads")
        paths.add_argument("-2", "--reverse", metavar="FASTQ", nargs='*', required=False, default=[],
            help="reverse paired-end reads in the same order")
        paths.add_argument("-s", "--single", metavar="FASTQ", nargs='*', required=False, default=[],
            help="single-end reads")
        paths.add_argument("-i", "--interleaved", metavar="FASTQ", nargs='*', required=False, default=[],
            help="interleaved reads")
        paths.add_argument("-b", "--background", metavar="FASTA", required=False,
            help="background")
        paths.add_argument("-o", "--output", metavar="PATH", required=True,
            help="path to output folder, will be created if non-existent")
        
        # "options" group
        parser.add_argument("-a", "--assemblers", nargs='*', required=False, default=ASSEMBLERS[:2],
            help=f"assemblers to use, pick any combination of {ASSEMBLERS}")
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
        # verify inputs
        #########################
        input_error = False
        def _error(message: str):
            nonlocal input_error
            parser.print_help()
            print(f"Invalid input: {message}")
            input_error = True

        output = Path(args.output).absolute()
        logs = output.joinpath("logs")
        if not output.exists(): os.makedirs(output)
        elif args.overwrite: shutil.rmtree(output)
        if not logs.exists(): os.makedirs(logs)

        _listify = lambda arg: [Path(p).absolute() for p in arg]
        reads_manifest = ReadsManifest(
            forward=_listify(args.forward),
            reverse=_listify(args.reverse),
            interleaved=_listify(args.interleaved),
            single=_listify(args.single),
        )
        reads_manifest.Verify(_error)
        reads_save = output.joinpath(ORIGINAL_READS)
        if not reads_save.parent.exists(): os.makedirs(reads_save.parent)
        reads_manifest.Save(reads_save)

        smk_args = ["--latency-wait 0"]
        for a in args.snakemake:
            if "=" in a:
                toks = a.split("=")
                pa = f"--{toks[0]} {'='.join(toks[1:])}"
            else:
                pa = f"--{a}"
            smk_args.append(pa)

        picked_assemblers = []
        _seen = set()
        for a in args.assemblers:
            if a in _seen: continue
            if a not in ASSEMBLERS:
                _error(f"[{a}] is not one of {ASSEMBLERS}")
                break
            picked_assemblers.append(a)
            _seen.add(a)

        if input_error: return
        smk_log = logs.joinpath("snakemake")
        link_log = "" if smk_log.exists() else f'ln -s {output.joinpath(".snakemake/log")} {logs.joinpath("snakemake")}'
        params = dict(
            src=MODULE_ROOT,
            log=logs,
            threads=args.threads,
            assemblers=','.join(picked_assemblers),
        )
        if args.background is not None:
            bg = Path(args.background).absolute()
        else:
            bg = output.joinpath(SKIP_FILTER)
            os.system(f"touch {bg}")
        params |= dict(background=bg)

        params_str = ' '.join(f"{k}={v}" for k, v in params.items())
        cmd = f"""\
            {link_log}
            snakemake -s {MODULE_ROOT.joinpath('main.smk')} -d {output}\
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
        mo.Procedure(args.args)

    def help(self, args=None):
        help = [
            f"{NAME} v{VERSION}",
            f"https://github.com/USER/{NAME}",
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
        cli.help # default
    )(cli, sys.argv[2:])

if __name__ == "__main__":
    main()
