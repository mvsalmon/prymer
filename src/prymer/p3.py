"""Primer design and related functions"""

import json
import logging
from pathlib import Path
import re
from typing import Any

import primer3

# regex for primer3 output parsing
P3_KEY = re.compile(r"PRIMER_(?P<kind>LEFT|RIGHT|PAIR|INTERNAL)_(?P<idx>\d+)(?:_(?P<field>.+))?")
# primer3 config file expected sections. Add sequence tags later?
P3_CONFIG_SECTIONS = {"global_tags"}

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

def default_product_size_range(target_len: int, global_tags: dict, slack_bp: int = 100) -> list[int]:
    """Calculate the product size range [min, max] based on input coordinates when a same chrom range is given
    Minimum: target_len + shortest primer on each side
    Maximum: Minimum + slack_bp per side (to avoid filling the whole template)."""
    min_primer = global_tags.get("PRIMER_MIN_SIZE", 18)
    lower = target_len + 2 * min_primer
    return [lower, lower + 2 * slack_bp]

def load_primer3_config(path: Path) -> dict[str, Any]:
    """Load primer3 configuration from file.
    Expected structure {"global_tags": {...}}.
    Parameters:
        path: Path to primer3 configuration file.
    Returns:
        global_tags: Dictionary of primer3 global tag:value pairs."""
    with open(path) as f:
        try:
            data = json.load(f)
        except json.decoder.JSONDecodeError as err:
            raise ValueError(f"{path} is not a valid JSON file.") from err

    if not isinstance(data, dict):
        raise ValueError(f"{path}: top level must be a valid object with 'global_tags' key.")
    unknown = set(data) - P3_CONFIG_SECTIONS
    if unknown:
        raise ValueError(f"{path}: unknown section(s): {sorted(unknown)};"
                         f"supported: {sorted(P3_CONFIG_SECTIONS)}.")
    global_tags = data.get("global_tags", {})
    if not isinstance(global_tags, dict):
        raise ValueError(f"{path}: 'global_tags' must be a dictionary object of TAG: value pairs.")

    if len(global_tags) == 0:
        raise ValueError(f"{path}: global_tags must not be empty. "
                         f"Add configurations or re-run without specifying --p3-params.")

    return global_tags

def parse_primer3_opts(p3_opts: list) -> dict[str, Any]:
    """Parse inputs for primer3 global and sequence tags and return a dictionary containing primer3 options with
    appropriate formatting: Ranges should be converted to a list, lists of ranges to lists of lists.
    each option should be in the format 'TAG=VALUE' Don't need to validate arguments here as primer3 will complain about
    invalid options.
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
        template: String. DNA sequence
        seq_tags: Dictionary of primer3 sequence tag:value pairs.
        global_tags: Dictionary of primer3 global tag:value pairs.
    Returns:
        Dictionary containing Primer3 output
    """
    logger.info("Designing primers...")
    seq_args = {**seq_tags, "SEQUENCE_TEMPLATE": template}
    return primer3.design_primers(seq_args=seq_args, global_args=global_tags)
