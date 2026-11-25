"""Functions for checking template region for SNPs"""


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

def snp_exclusion_regions(records: list[dict], win_start: int, win_end: int) -> list[list[int]]:
    """Convert UCSC variant records to merged template-relative [start, length] regions.
    Output passed to primer3 SEQUENCE_EXCLUDED_REGIONS"""
    regions = []
    for rec in records:
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