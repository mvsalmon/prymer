"""Functions handling commandline interface, argparse etc."""

import argparse
import logging
from pathlib import Path

from prymer.coordinates import parse_coordinates
from prymer.exceptions import UCSCError
from prymer.models import Coordinate, Job, Settings
from prymer.output import validate_output_path, write_outputs
from prymer.p3 import parse_primer3_opts, load_primer3_config
from prymer.pipeline import run_job
from prymer.snps import load_maf_manifest
from prymer.ucsc import fetch_sequence

# Argument defaults
DEFAULTS = {
    "template_sequence_length": 500,
    "reference_genome": "hg38",
    "output_path": "./",
    "breakpoint_flanking_bp": 50,
    "start_primer_position": "5",
    "end_primer_position": "3",
    "snp_freq_cutoff": 0.01,
    "maf_source": "Any"
}

logger = logging.getLogger(__name__)

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="prymer",
        description="Design PCR primers for given genomic coordinates"
    )
    parser.add_argument(
        "-s", "--start-coordinate", type=coordinate_arg, required=True,
        help="""Genomic coordinate e.g "chr1:23,456" (quoted if it has commas). 
                With --sv-breakpoint, this is the 5' breakpoint""",
   )
    parser.add_argument(
        "-e", "--end-coordinate", type=coordinate_arg,
        help="Optional second genomic coordinate. Specifies the end of the PCR target range."
             "When --sv-breakpoint is used, this is mandatory and specifies the 3' breakpoint location."
    )
    parser.add_argument(
        "-l", "--template-sequence-length", type=int,
        default=DEFAULTS["template_sequence_length"],
        help="Length of template sequence (in bp) around the target sequence in which primers will be designed." 
             "For SV, this is per-side of the breakpoint. Default: %(default)s"
    )
    parser.add_argument(
        "-r", "--reference-genome", default=DEFAULTS["reference_genome"],
        help="UCSC genome name. Default: %(default)s."
    )
    parser.add_argument(
        "-f", "--sv-breakpoint", action="store_true",
        help="Design primers spanning an SV breakpoint. i.e. --start-cordinate and --end-cordinate are two regions that are normally non-contiguous",
    )
    parser.add_argument(
        "-o", "--output-path", type=Path, default=DEFAULTS["output_path"],
        help="Output filepath. Defaults: current directory",
    )
    parser.add_argument(
        "-n", "--output-name",
        help="Output file name (no extension). Default derived from coordinates"
    )
    parser.add_argument(
        "--start-primer-position", choices=["5", "3"],
        default=DEFAULTS["start_primer_position"],
        help="""Specify the position of the primer relative to the breakpoint (5' or 3') for the 
                        --start_coordinate. Choose one of 5 or 3. Default %(default)s. 
                        Note: do not include the trailing apostrophe."""
    )
    parser.add_argument(
        "--end-primer-position", choices=["5", "3"],
        default=DEFAULTS["end_primer_position"],
        help="As --start_primer_position, for --end_coordinate. Default %(default)s.",
    )
    parser.add_argument(
        "--reverse-complement", choices=["start", "end", "both"],
        help="Reverse complement the template sequence(s)."
    )
    parser.add_argument(
        "--p3-global-tags", nargs="*", default=[], metavar="TAG=VALUE",
        help="""Specify additional primer3 global tags as a space separated list in the form "OPTION1_NAME=<value(s)>" "OPTION2_NAME=<value(s)>". 
                The version of primer3 used is 2.6.1. See https://htmlpreview.github.io/?https://github.com/primer3-org/primer3/blob/v2.6.1/src/primer3_manual.htm
                for details of all options available.""",
    )
    parser.add_argument(
        "--p3-sequence-tags", nargs="*", default=[], metavar="TAG=VALUE",
        help="""Specify additional primer3 sequence tags as a space separated list in the form "OPTION1_NAME=<value(s)>" "OPTION2_NAME=<value(s)>".
         The version of primer3 used is 2.6.1. See https://htmlpreview.github.io/?https://github.com/primer3-org/primer3/blob/v2.6.1/src/primer3_manual.htm
        for details of all options available.""",
    )
    parser.add_argument(
        "--snp-track", type=str,
        help="UCSC dbSNP track, e.g. dbSnp155Common. If set, primers will not be placed over variants "
             "present in this track. Can be used with --snp-freq-cutoff."
    )
    parser.add_argument(
        "--snp-freq-cutoff", type=float, default=None,
        help="If --snp-track is specified, SNPs in variant track with a frequency below this level will be ignored, and may be included"
             "in the primer binding region. Must be used in conjunction with --snp-track (default: 0.01)."
    )
    parser.add_argument(
        "--maf-source", type=str, nargs="*", default=None,
        help="Define 1 or more sources for the dbSNP MAF value, as specified with --snp-maf-cutoff. If no source is given, "
             "a SNP will be excluded from primer regions if the MAF from any source is above the cutoff level (default 0.01)."
             "Must be used in conjunction with --snp-track"
    )
    parser.add_argument(
        "--p3-config", type=Path, metavar="FILE",
        help="Path to JSON file containing preset primer3 global tag configuration options."
             "Options given with --p3-global-tags will override those in the file."
             " See p3_example_params.json for the required file structure."
    )
    return parser

def coordinate_arg(text: str) -> Coordinate:
    """argparse 'type=' adapter: turns ValueError into clean usage error."""
    try:
        return parse_coordinates(text)
    except ValueError as err:
        raise argparse.ArgumentTypeError(str(err)) from err

def default_job_name(start: Coordinate, end: Coordinate | None) -> str:
    """Generate a default job name from input coordinates."""
    name = f"{start.chrom}_{start.pos}"
    return name if end is None else f"{name}_{end.chrom}_{end.pos}"

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse and check input args."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.sv_breakpoint and args.end_coordinate is None:
        parser.error("--end-coordinate is required with --sv-brakpoint")
    if args.reverse_complement and not args.sv_breakpoint:
        parser.error("--reverse_complement only allowed with --sv-breakpoint.")
    if args.template_sequence_length < 1:
        parser.error("--template_sequence_length must be a positive integer.")

    # SNP checks
    if args.snp_freq_cutoff is not None and args.snp_track is None:
        parser.error("--snp-freq-cutoff only allowed with --snp-track.")
    if args.maf_source is not None and args.snp_track is None:
        parser.error("--maf-source only allowed with --snp-track.")

    # Default: no snp check
    args.maf_manifest = None

    if args.snp_track:
        if args.snp_freq_cutoff is None:
            args.snp_freq_cutoff = DEFAULTS["snp_freq_cutoff"]

    if args.maf_source:
        manifest = load_maf_manifest()
        unknown = [s for s in args.maf_source if s not in manifest]
        if unknown:
            parser.error(f"Unknown MAF source(s): {', '.join(unknown)}. "
                         "Run with --list-maf-sources for valid options.")
        args.maf_manifest = manifest
    else:
        args.maf_source = None

    # Check for primer3 params file
    args.p3_file_global = {}
    if args.p3_config:
        try:
            args.p3_file_global = load_primer3_config(args.p3_config)
        except (OSError, ValueError) as err:
            parser.error(str(err))

    return args

def settings_from_args(args: argparse.Namespace) -> Settings:
    """Generate Settings object from arguments"""
    return Settings(
        ref_genome=args.reference_genome,
        output_dir=validate_output_path(args.output_path),
        p3_global_tags={**args.p3_file_global, **parse_primer3_opts(args.p3_global_tags)},
        snp_track=args.snp_track,
        snp_freq_cutoff=args.snp_freq_cutoff,
        maf_source=args.maf_source,
        maf_manifest=tuple(args.maf_manifest) if args.maf_manifest else None,
    )

def jobs_from_args(args: argparse.Namespace) -> list[Job]:
    """Generate Job objects from arguments.
    Returns a list of Job objects"""
    return [
        Job(
            name=args.output_name or default_job_name(args.start_coordinate, args.end_coordinate),
            start=args.start_coordinate,
            end=args.end_coordinate,
            sv=args.sv_breakpoint,
            start_primer_position=args.start_primer_position,
            end_primer_position=args.end_primer_position,
            rev_comp=args.reverse_complement,
            seq_len=args.template_sequence_length,
            p3_seq_tags=parse_primer3_opts(args.p3_sequence_tags),
        )
    ]

def main(argv: list[str] | None = None, fetch=fetch_sequence) -> int:
   logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
   args = parse_args(argv)
   try:
       settings = settings_from_args(args)
       jobs = jobs_from_args(args)
       results = [run_job(job, settings, fetch=fetch) for job in jobs]
       write_outputs(results, settings)
   except (ValueError, UCSCError, OSError) as err:
       logger.error("%s", err)
       return 1
   return 0
