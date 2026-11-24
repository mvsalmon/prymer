"""Run primer design jobs"""

import logging

from prymer.coordinates import (
    breakpoint_window,
    single_chrom_range_window,
    single_locus_window,
)
from prymer.models import Job, Result, Settings
from prymer.p3 import default_product_size_range, parse_primer3_output, primer_design
from prymer.sequence import build_breakpoint
from prymer.snps import add_excluded_regions, snp_exclusion_regions
from prymer.ucsc import build_ucsc_track_url, build_ucsc_url, fetch_sequence, fetch_snps

logger = logging.getLogger(__name__)

def run_job(job: Job, settings: Settings, fetch=fetch_sequence) -> Result:
    """Runs primer design job.
    Parameters:
        job: Job object
        settings: Settings object
        fetch: Dependency injection of fetch_sequence function.
    Returns:
        Result object containing parsed primer3 output"""

    global_tags = settings.p3_global_tags

    # Primers for single locus (i.e. non-fusion/SV)
    if not job.sv and job.end is None:
        template_win = single_locus_window(job.start, job.seq_len)
        template = fetch(build_ucsc_url(settings.ref_genome, *template_win))

        # Target is 1bp in middle of the template
        seq_target_pos = len(template) // 2
        seq_tags = {
            "SEQUENCE_TARGET": [seq_target_pos, 1],
            **job.p3_seq_tags,
        }
    # Primers where start and end specify the template range
    elif not job.sv and job.end is not None:
        template_win = single_chrom_range_window(job.start, job.end, job.seq_len)
        template = fetch(build_ucsc_url(settings.ref_genome, *template_win))
        # Target is the range between start and end pos which begins after job.seq_len bp
        seq_target_len = job.end.pos - job.start.pos
        seq_tags = {
            "SEQUENCE_TARGET": [job.seq_len, seq_target_len],
            **job.p3_seq_tags,
        }
        # Calculate Primer product size range based on target size
        product_size_range = default_product_size_range(seq_target_len, settings.p3_global_tags)
        global_tags = {
            "PRIMER_PRODUCT_SIZE_RANGE": product_size_range,
            **settings.p3_global_tags,
        }
        # Check range size is valid
        rng = global_tags["PRIMER_PRODUCT_SIZE_RANGE"]
        ranges = rng if isinstance(rng[0], list) else [rng]
        min_possible = seq_target_len + 2 * global_tags.get("PRIMER_MIN_SIZE", 18)
        if max(r[1] for r in ranges) < min_possible:
            raise ValueError(f"PRIMER_PRODUCT_SIZE_RANGE max is below the minimum {min_possible}bp needed "
                             f"to span a {seq_target_len}bp target.")
    else:
        assert job.end is not None

        start_win = breakpoint_window(job.start, job.start_primer_position, job.seq_len)
        end_win = breakpoint_window(job.end, job.end_primer_position, job.seq_len)

        start_seq = fetch(build_ucsc_url(settings.ref_genome, *start_win))
        end_seq = fetch(build_ucsc_url(settings.ref_genome, *end_win))

        template = build_breakpoint(start_seq, end_seq, job.rev_comp)

        junction = len(start_seq)
        seq_tags = {
            "SEQUENCE_TARGET": [max(0, junction - job.flank_bp), 2 * job.flank_bp],
            **job.p3_seq_tags,
        }
    # SNP checks if asked for
    if settings.snp_track:
        if job.sv:
            logger.warning("%s: SNP filtering not supported for SV jobs; skipping. "
                           "Check primer binding regions for SNPS manually if needed.", job.name)
        else:
            _, win_start, win_end = template_win
            url = build_ucsc_track_url(settings.ref_genome, settings.snp_track, *template_win)
            records = fetch_snps(url, settings.snp_track)
            seq_tags = add_excluded_regions(seq_tags, snp_exclusion_regions(records, win_start, win_end))
            logger.info("%s: excluding %d SNP region(s)", job.name, len(seq_tags.get("SEQUENCE_EXCLUDED_REGION", [])))

    raw = primer_design(template, seq_tags, global_tags)
    return Result(job, parse_primer3_output(raw), raw["PRIMER_PAIR_EXPLAIN"])