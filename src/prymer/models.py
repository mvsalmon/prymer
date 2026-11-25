"""Shared data structures for Prymer"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

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
        if self.end is not None and not self.sv and self.end.pos < self.start.pos:
            raise ValueError(f"Start position {self.start.pos} falls after end position {self.end.pos}; "
                             f"start must be before end position for same chromosome jobs."  )
        if not self.sv and self.end is not None and self.start.chrom != self.end.chrom:
            raise ValueError(f"{self.start.chrom} and {self.end.chrom} differ; use --sv-breakpoint "
                             "or check coordinates.")



@dataclass(frozen=True)
class Settings:
    """Dataclass to hold primer design settings."""
    ref_genome: str
    output_dir: Path
    p3_global_tags: dict
    snp_track: str | None = None

@dataclass(frozen=True)
class Result:
    """Dataclass to hold results of primer design job."""
    job: Job
    pairs: dict[str, dict]
    explain: str
