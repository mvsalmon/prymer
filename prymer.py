# TODO take file with coordinates as input for multi job runs
# TODO handle a design job where the input is a single range, not just a single nucleotide
# TODO add PRIMER_EXPLAIN info to output file?

import argparse
import logging
import os
import re
import json
from pathlib import Path
from typing import Optional, Dict, NamedTuple, Tuple, Any
from dataclasses import dataclass, field
import requests
import pandas as pd
import primer3
from primer3 import p3helpers

class Coordinate(NamedTuple):
    chrom: str
    pos: int

@dataclass(frozen=True)
class Job:
    """Dataclass to represent a single primer design job."""
    name: str
    start: Coordinate
    end: Coordinate | None = None
    sv: bool = False
    start_primer_position: str = "5"
    end_primer_position: str = "3"
    rev_comp: str | None = None
    seq_len: int = 500
    flank_bp: int = 50
    p3_seq_tags: dict = field(default_factory=dict)

    def __post_init__(self):
        # Validation of various Job params
        if self.sv and self.end is None:
            raise ValueError("SV design jobs need an end coordinate.")
        if self.start.pos - self.seq_len < 0:
            raise ValueError(
                f"Template sequence start = {self.start.pos - self.seq_len}. Must be greater than 0; "
                "try a shorter --template-sequence-length.")
        if not self.sv and self.end is not None and self.start.chrom != self.end.chrom:
            raise ValueError(f"{self.start.chrom} and {self.end.chrom} differ; use --sv-breakpoint "
                             "or check coordinates.")


@dataclass(frozen=True)
class Settings:
    """Dataclass to hold primer design settings."""
    ref_genome: str
    output_dir: Path
    p3_global_tags: dict

@dataclass(frozen=True)
class Result:
    """Dataclass to hold results of primer design job."""
    job: Job
    pairs: dict[str, dict]
    explain: str

class UCSCError(Exception):
    """Raised when UCSC request fails."""

# Argument defaults
DEFAULTS = {
    "template_sequence_length": 500,
    "reference_genome": "hg38",
    "output_path": "./",
    "breakpoint_flanking_bp": 50,
    "start_primer_position": "5",
    "end_primer_position": "3"
}

# regex for primer3 output parsing
P3_KEY = re.compile(r"PRIMER_(?P<kind>LEFT|RIGHT|PAIR|INTERNAL)_(?P<idx>\d+)(?:_(?P<field>.+))?")

logger = logging.getLogger(__name__)

def _coerce(token: str):
    """Attempt to coerce a string to an int or float if possible, else leave as str."""
    for caster in (int, float):
        try:
            return caster(token)
        except ValueError:
            pass
    return token

def _parse_item(item: str):
    """Parse one item from a primer3 option:
    '20' -> 20, '50,400' -> [50, 400], '100-300' -> [100, 300]."""
    if "," in item:
        return [_coerce(x) for x in item.split(",")]
    # size ranges
    if re.fullmatch(r"\d+-\d+", item):
        return [int(x) for x in item.split("-")]
    # single value
    return _coerce(item)

def validate_output_path(path_str: str) -> Path:
    """Validate the output directory exists and is writable.
    Parameters:
        path_str: String containing user-specified output directory path
    Returns:
          Output directory path.
    Raises:
        FileNotFoundError: If output directory does not exist.
        PermissionError: If output directory is not writable.
        """
    path = Path(path_str)
    if not path.exists():
        raise FileNotFoundError(f"Output path does not exist: {path}")
    if not os.access(path, os.W_OK):
        raise PermissionError(f"Output path is not writable: {path}")
    return path

def parse_coordinates(coordinate: str) -> Coordinate[str, int]:
    """Parse and validate coordinates (e.g. 'chr1:123456').
    Parameters:
        coordinate: String containing coordinate to be validated.
    Returns:
        Tuple containing the chromosome (str) and position (int), e.g. ['chr1', 123456].
    Raises:
        ValueError: If coordinate is invalid.
    """
    clean_coordinate = coordinate.replace(",", "")

    if not re.fullmatch(r'chr([1-9]|1[0-9]|2[0-2]|[XYM]):\d+', clean_coordinate):
        raise ValueError(f"Invalid coordinates '{coordinate}'. Must be in the format 'chr1:123456'.")

    chrom, pos = clean_coordinate.split(":")
    return Coordinate(chrom, int(pos))

def coordinate_arg(text: str) -> Coordinate:
    """argparse 'type=' adapter: turns ValueError into clean usage error."""
    try:
        return parse_coordinates(text)
    except ValueError as err:
        raise argparse.ArgumentTypeError(str(err)) from err

def parse_primer3_opts(p3_opts: list) -> Dict[str, Any]:
    """Parse inputs for primer3 global and sequence tags and return a dictionary containing primer3 options with
    appropriate formatting: Ranges should be converted to a list, lists of ranges to lists of lists.
    each option should be in the format 'TAG=VALUE'.
    Parameters:
        p3_opts: list of p3 tag:value pairs.
    Returns:
        Dictionary of primer3 input options.
    """
    p3_options = {}
    for option in p3_opts:

        if "=" not in option:
            raise ValueError(f"Missing '=' in primer3 option '{option}'.")

        tag, _, value = option.partition("=")
        # parse values with multiple options, split on non-alphanumeric chars.
        items = [x for x in re.split(r"[;\s]+", value.strip()) if x]
        parsed = [_parse_item(x) for x in items]

        p3_options[tag.strip()] = parsed[0] if len(parsed) == 1 else parsed

    return p3_options

    # p3helpers only included in v1.2.0+
def reverse_complement(sequence: str) -> str:
    """Wrapper around primer3.p3helpers.reverse_complement()."""
    return primer3.p3helpers.reverse_complement(sequence)

def parse_primer3_output(raw: dict) -> dict[str, dict]:
    """Parse raw primer3 output. Ignores PRIMER_PAIR_EXPLAIN which can be accessed directly.
    Uses P3_KEY regex"""
    pairs: dict[str, dict] = {}
    for key, value in raw.items():
        m = P3_KEY.fullmatch(key)
        if m is None:
            continue
        name = f"PRIMER_{m['kind']}" + (f"_{m['field']}" if m["field"] else "")
        pairs.setdefault(f"PAIR_{m['idx']}", {})[name] = value
    return pairs

def primer_design(template: str, seq_tags: dict, global_tags: dict) -> dict:
    """Design PCR primers using primer3 with default options
    Parameters:
        primer3_options:
    Returns:
        Dictionary containing Primer3 output
    """
    logger.info("Designing primers...")
    seq_args = {**seq_tags, "SEQUENCE_TEMPLATE": template}
    return primer3.design_primers(seq_args=seq_args, global_args=global_tags)

def build_UCSC_url(genome: str, chrom: str, start: int, end: int) -> str:
    """Construct URL for UCSC API request."""
    return(f"https://api.genome.ucsc.edu/getData/sequence?"
           f"genome={genome};chrom={chrom};start={start};end={end}")

def fetch_sequence(url: str, timeout: int = 10) -> str:
    """Fetch sequence from UCSC API. Returns string containing target sequence."""
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
    except requests.RequestException as err:
        raise UCSCError(f"UCSC request failed: {err}\n {json.loads(err.response.text)["error"]}") from err
    return response.json()['dna'].upper()

def breakpoint_window(coordinate: Tuple[str, int], position: str, length: int) -> Tuple[str, int, int]:
    """Returns (chrom, start, end) for breakpoint template sequences.
    Sequence ends at the breakpoint if the primer position is 5' (default for start pos).
    Sequence starts at the breakpoint if the primer position is 3' (default for end pos).
    """
    chrom, pos = coordinate
    if position == "5":
        return chrom, pos - length, pos
    if position == "3":
        return chrom, pos, pos + length
    raise ValueError(f"Invalid breakpoint position: {position}")

def single_locus_window(coordinate: Tuple[str, int], length: int) -> Tuple[str, int, int]:
    """Returns (chrom, start, end) window for single locus templates.
    e.g. for a single coordinate 'chr1:1000' and length 500, we return ('chr1', 500, 1500)"""
    chrom, pos = coordinate
    return chrom, pos - length, pos + length

def single_chrom_range_window(start_coordinate: Coordinate, end_coordinate: Coordinate,
                              length: int) -> Tuple[str, int, int]:
    """Returns (chrom, start, end) window for target range on single chrom."""
    chrom, start_pos = start_coordinate
    end_pos = end_coordinate.pos
    return chrom, start_pos - length, end_pos + length

def build_breakpoint(start_seq: str, end_seq: str, rev_comp: str | None = None) -> str:
    """Join breakpoint sequences for fusion primers.
    Parameters:
        start_seq: String. Sequence of start (first) breakpoint.
        end_seq: String. Sequence of the end (second) breakpoint.
        rev_comp: String. Controls reverse complementation. Can be 'start', 'end', 'both', or None/empty.

    Returns:
        String of the full breakpoint sequence for primer design.
    """
    if rev_comp in ("start", "both"):
        start_seq = reverse_complement(start_seq)
    if rev_comp in ("end", "both"):
        end_seq = reverse_complement(end_seq)
    return start_seq + end_seq

def run_job(job: Job, settings: Settings, fetch=fetch_sequence) -> Result:
    """Runs primer design job.
    Parameters:
        job: Job object
        settings: Settings object
        fetch: Dependency injection of fetch_sequence function.
    Returns:
        Result object containing parsed primer3 output"""
    # Primers for single locus (i.e. non-fusion/SV)
    if not job.sv and job.end is None:
        template_win = single_locus_window(job.start, job.seq_len)
        template = fetch(build_UCSC_url(settings.ref_genome, *template_win))

        # Target is 2bp in middle of the template
        seq_target_pos = len(template) // 2
        seq_tags = {
            "SEQUENCE_TARGET": [seq_target_pos, 1],
            **job.p3_seq_tags,
        }
    # Primers where start and end specify the template range
    elif not job.sv and job.end is not None:
        template_win = single_chrom_range_window(job.start, job.end, job.seq_len)
        template = fetch(build_UCSC_url(settings.ref_genome, *template_win))

        # Target is the range between start and end pos which begins after job.seq_len bp
        seq_target_len = job.end.pos - job.start.pos
        seq_tags = {
            "SEQUENCE_TARGET": [job.seq_len, seq_target_len],
            **job.p3_seq_tags,
        }
    else:
        assert job.end is not None

        start_win = breakpoint_window(job.start, job.start_primer_position, job.seq_len)
        end_win = breakpoint_window(job.end, job.end_primer_position, job.seq_len)

        start_seq = fetch(build_UCSC_url(settings.ref_genome, *start_win))
        end_seq = fetch(build_UCSC_url(settings.ref_genome, *end_win))

        template = build_breakpoint(start_seq, end_seq, job.rev_comp)

        junction = len(start_seq)
        seq_tags = {
            "SEQUENCE_TARGET": [max(0, junction - job.flank_bp), 2 * job.flank_bp],
            **job.p3_seq_tags,
        }

    raw = primer_design(template, seq_tags, settings.p3_global_tags)
    return Result(job, parse_primer3_output(raw), raw["PRIMER_PAIR_EXPLAIN"])

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
        help=f"""Specify the position of the primer relative to the breakpoint (5' or 3') for the 
                        --start_coordinate. Choose one of 5 or 3. Default %(default)s. 
                        Note: do not include the trailing apostrophe."""
    )
    parser.add_argument(
        "--end-primer-position", choices=["5", "3"],
        default=DEFAULTS["end_primer_position"],
        help=f"As --start_primer_position, for --end_coordinate. Default %(default)s.",
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
    return parser

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
    if args.start_coordinate and args.end_coordinate and args.sv_breakpoint is None:
        # Catch this for now, want to implement a range input later
        parser.error("Using --start-coordinate and --end-coordinate together requires --sv-breakpoint.")
    return args

def default_job_name(start: Coordinate, end: Coordinate | None) -> str:
    """Generate a default job name from input coordinates."""
    name = f"{start.chrom}_{start.pos}"
    return name if end is None else f"{name}_{end.chrom}_{end.pos}"

def settings_from_args(args: argparse.Namespace) -> Settings:
    """Generate Settings object from arguments"""
    return Settings(
        ref_genome=args.reference_genome,
        output_dir=validate_output_path(args.output_path),
        p3_global_tags=parse_primer3_opts(args.p3_global_tags)
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

def write_output(result: Result, output_dir: Path) -> Path | None:
    """Write results of one Job to disk. Returns the path or None if no primers."""
    logger.info("%s: %s", result.job.name, result.explain)
    if not result.pairs:
        logger.warning("%s: Primer3 found no valid primers", result.job.name)
        return None
    path = output_dir / f"{result.job.name}.csv"
    pd.DataFrame.from_dict(result.pairs).to_csv(path)
    return path

def write_outputs(results: list[Result], settings: Settings) -> list[Path]:
    paths = (write_output(r, settings.output_dir) for r in results)
    return [p for p in paths if p is not None]


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


if __name__ == "__main__":
    raise SystemExit(main())
