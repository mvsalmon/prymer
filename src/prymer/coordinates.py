"""Genomic coordinate handling"""

import re

from prymer.models import Coordinate


def parse_coordinates(coordinate: str) -> Coordinate:
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

def breakpoint_window(coordinate: Coordinate, position: str, length: int) -> tuple[str, int, int]:
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

def single_locus_window(coordinate: Coordinate, length: int) -> tuple[str, int, int]:
    """Returns genomic (chrom, start, end) window for single locus templates.
    e.g. for a single coordinate 'chr1:1000' and length 500, we return ('chr1', 500, 1500)"""
    chrom, pos = coordinate
    return chrom, pos - length, pos + length

def single_chrom_range_window(start_coordinate: Coordinate, end_coordinate: Coordinate,
                              length: int) -> tuple[str, int, int]:
    """Returns (chrom, start, end) window for target range on single chrom."""
    chrom, start_pos = start_coordinate
    end_pos = end_coordinate.pos
    return chrom, start_pos - length, end_pos + length
