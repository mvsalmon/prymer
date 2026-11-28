"""Functions related to checking template region for SNPs"""
import logging
from importlib.resources import files
from pathlib import Path

from prymer.models import Settings

logger = logging.getLogger(__name__)

def _maf_source_index(maf_source: str, maf_manifest: tuple[str, ...]) -> int:
    """Returns the index of the MAF source in MAF manifest"""
    return maf_manifest.index(maf_source)

def load_maf_manifest(path: Path | None = None) -> list[str]:
    """Load contents of MAF manifest file. This contains the sources of MAF values contained in the UCSC
    dbSNP155 track, in the order they appear.
    Parameters:
        path: Path to MAF manifest file
    Returns:
        List containing the MAF value sources."""
    if path is None:
        text = (files("prymer") / "data" / "MAF.txt").read_text(encoding="utf-8")
    else:
        text = Path(path).read_text(encoding="utf-8")
    return [
        line.rstrip('\n')
        for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")]


def check_maf_threshold(var_record: dict,
              maf_threshold: float,
              maf_manifest: tuple[str, ...] | None = None,
              maf_source: list[str] | None = None,
              ) -> bool:
    """Check MAF from specified source against threshold cutoff. Returns True if the MAF is above the threshold
    and should be excluded from primer design. Otherwise, return False."""
    mafs = var_record.get("minorAlleleFreq", [])
    # Decide what to do for missing maf data. Ignore position for now - low evidence for true SNP.
    if not mafs:
        logger.warning(f"No MAF data found for {var_record['name']}. Variant class is {var_record['class']}. Position will be included.")
        return False
    # MAF values are comma separated in a single string and have a trailing comma - split and convert to tuple of floats
    mafs = tuple(float(x) for x in mafs.split(',') if x.strip())
    if maf_source is None:
        # True if any value in mafs > threshold when no specific source given.
        return any(float(maf) > maf_threshold for maf in mafs)

    # True if maf > threshold in any source
    # Update this to show which source and MAF were above threshold at some point.
    maf_indices = [_maf_source_index(source, maf_manifest) for source in maf_source]
    if any(float(mafs[i]) > maf_threshold for i in maf_indices):
        logger.info(f"MAF for {var_record['name']} above threshold. Position will be excluded.")
        return True
    else:
        return False

def merge_regions(regions: list[list[int]]) -> list[list[int]]:
    """Merge overlapping or adjacent [start, length] variant regions.
    Sort by start, then if then next region begins at or before the end of the
    current one, extend the current one; otherwise start a new one"""
    merged: list[list[int]] = []
    for start, length in sorted(regions):
        end = start + length
        if merged and start <= merged[-1][0] + merged[-1][1]:
            # prev is just a reference to merged, so mutating prev also mutates merged
            prev = merged[-1]
            prev[1] = max(prev[0] + prev[1], end) - prev[0]
        else: merged.append([start, length])
    return merged

def snp_exclusion_regions(records: list[dict], win_start: int, win_end: int, settings: Settings) -> list[list[int]]:
    """Convert UCSC variant records to merged template-relative [start, length] regions.
    Output passed to primer3 SEQUENCE_EXCLUDED_REGIONS"""
    regions = []
    for rec in records:
        # Ignore SNP if MAF below threshold
        if not check_maf_threshold(rec, settings.maf_threshold, settings.maf_manifest, settings.maf_source):
            continue
        snp_start, snp_end = rec["chromStart"], rec["chromEnd"]
        # handle 0-length insertions vars with flanking bp
        if snp_start == snp_end:
            snp_start, snp_end = snp_start - 1, snp_end + 1
        # clip variants that span template window edges
        snp_start, snp_end = max(snp_start, win_start), min(snp_end, win_end)
        if snp_end > snp_start:
            regions.append([snp_start - win_start, snp_end - snp_start])
    return merge_regions(regions)

def add_excluded_regions(seq_tags: dict, regions: list[list[int]]) -> dict:
    """Combine SNP regions with user-specified excluded regions"""
    existing = seq_tags.get("SEQUENCE_EXCLUDED_REGION", [])
    if existing and not isinstance(existing[0], list):
        existing = [existing]
    merged = merge_regions([*existing, *regions])
    return {**seq_tags, "SEQUENCE_EXCLUDED_REGION": merged} if merged else seq_tags

