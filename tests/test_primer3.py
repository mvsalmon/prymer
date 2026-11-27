import json
import re
from pathlib import Path

import primer3
import pytest

from prymer.cli import parse_args, settings_from_args
from prymer.models import Coordinate, Job, Settings
from prymer.p3 import (
    load_primer3_config,
    parse_primer3_opts,
    parse_primer3_output,
    primer_design,
)
from prymer.pipeline import run_job


def test_comma_range():
    result = parse_primer3_opts(["PRIMER_PRODUCT_SIZE_RANGE=100-300 300-500"])
    assert isinstance(result["PRIMER_PRODUCT_SIZE_RANGE"], list)

@pytest.mark.parametrize(
    "options, expected",
    [
        (["PRIMER_OPT_TM=60.5"],  {"PRIMER_OPT_TM": 60.5}),
        (["PRIMER_MIN_THREE_PRIME_DISTANCE=-1"],  {"PRIMER_MIN_THREE_PRIME_DISTANCE": -1}),
        (["SEQUENCE_INCLUDED_REGION=50,400"],  {"SEQUENCE_INCLUDED_REGION": [50, 400]}),
        (["SEQUENCE_TARGET=50,2 100,2"],   {"SEQUENCE_TARGET": [[50, 2], [100, 2]]}),
        (["PRIMER_PRODUCT_SIZE_RANGE=100,300 400,500"],{"PRIMER_PRODUCT_SIZE_RANGE": [[100, 300], [400, 500]]}),
        (["SEQUENCE_ID=sample.1"], {"SEQUENCE_ID": 'sample.1'}),
        (["PRIMER_OPT_TM=60.5", "SEQUENCE_TARGET=50,2 100,2"], {"PRIMER_OPT_TM": 60.5, "SEQUENCE_TARGET": [[50, 2], [100,2]]})
    ],
)

def test_p3_option_parsing(options, expected):
    assert parse_primer3_opts(options) == expected

# Output parser tests
def test_groups_keys_by_pair_index():
    raw = {
        "PRIMER_PAIR_0_PENALTY": 0.1,
        "PRIMER_LEFT_0_SEQUENCE": "AAAA",
        "PRIMER_PAIR_1_PENALTY": 0.2,
        "PRIMER_LEFT_1_SEQUENCE": "CCCC",
    }
    assert parse_primer3_output(raw) == {
        "PAIR_0": {"PRIMER_PAIR_PENALTY": 0.1, "PRIMER_LEFT_SEQUENCE": "AAAA"},
        "PAIR_1": {"PRIMER_PAIR_PENALTY": 0.2, "PRIMER_LEFT_SEQUENCE": "CCCC"},
    }

def test_key_without_field_keeps_base_name():
    # PRIMER_LEFT_0 holds [start, length], with no trailing field name
    assert parse_primer3_output({"PRIMER_LEFT_0": [237, 20]}) == {
        "PAIR_0": {"PRIMER_LEFT": [237, 20]}
    }


def test_multi_digit_pair_index():
    # catches a regex using \d instead of \d+
    assert list(parse_primer3_output({"PRIMER_PAIR_10_PENALTY": 0.5})) == ["PAIR_10"]

@pytest.mark.parametrize("kind", ["LEFT", "RIGHT", "PAIR", "INTERNAL"])
def test_all_kinds_recognised(kind):
    result = parse_primer3_output({f"PRIMER_{kind}_0_PENALTY": 1.0})
    assert result == {"PAIR_0": {f"PRIMER_{kind}_PENALTY": 1.0}}

def test_empty_input_returns_empty_dict():
    assert parse_primer3_output({}) == {}


def test_no_primers_found_returns_empty_dict():
    # what primer3 returns when it can't design anything
    raw = {
        "PRIMER_PAIR_EXPLAIN": "considered 0, ok 0",
        "PRIMER_LEFT_NUM_RETURNED": 0,
        "PRIMER_PAIR_NUM_RETURNED": 0,
    }
    assert parse_primer3_output(raw) == {}

def test_one_entry_per_returned_pair(p3_raw_output):
    result = parse_primer3_output(p3_raw_output)
    n = p3_raw_output["PRIMER_PAIR_NUM_RETURNED"]
    assert set(result) == {f"PAIR_{i}" for i in range(n)}


def test_every_pair_has_identical_fields(p3_raw_output):
    # ragged pairs would become NaN cells in the output CSV
    result = parse_primer3_output(p3_raw_output)
    key_sets = [set(pair) for pair in result.values()]
    assert all(keys == key_sets[0] for keys in key_sets)


def test_expected_fields_present(p3_raw_output):
    pair = parse_primer3_output(p3_raw_output)["PAIR_0"]
    expected = {
        "PRIMER_LEFT_SEQUENCE", "PRIMER_RIGHT_SEQUENCE",
        "PRIMER_LEFT", "PRIMER_RIGHT",
        "PRIMER_LEFT_TM", "PRIMER_RIGHT_TM",
        "PRIMER_PAIR_PENALTY", "PRIMER_PAIR_PRODUCT_SIZE",
    }
    assert expected <= pair.keys()          # subset: extra keys are fine


def test_no_pair_index_left_in_field_names(p3_raw_output):
    result = parse_primer3_output(p3_raw_output)
    for pair in result.values():
        assert not any(re.search(r"\d", name) for name in pair)


def test_known_values_for_pair_0(p3_raw_output):
    # "golden" values copied from tests/data/p3_raw_out.json
    pair = parse_primer3_output(p3_raw_output)["PAIR_0"]
    assert pair["PRIMER_LEFT_SEQUENCE"] == "tcagcaggattttgacggct"
    assert pair["PRIMER_RIGHT_SEQUENCE"] == "ATGGactctgccagatccct"
    assert pair["PRIMER_LEFT"] == [237, 20]
    assert pair["PRIMER_PAIR_PRODUCT_SIZE"] == 172

def test_design_primers_sends_template(monkeypatch):
    captured = {}
    def fake(seq_args, global_args):
        captured.update(seq_args)
        return {"PRIMER_PAIR_EXPLAIN": "ok"}
    monkeypatch.setattr(primer3, "design_primers", fake)

    primer_design("AGCT", {"SEQUENCE_ID": "a"}, {})
    assert captured["SEQUENCE_TEMPLATE"] == "AGCT"

def test_range_target_sets_product_size(monkeypatch):
    captured = {}
    def fake(seq_args, global_args):
        captured.update(global_args)
        return {"PRIMER_PAIR_EXPLAIN": "ok"}
    monkeypatch.setattr(primer3, "design_primers", fake)

    job = Job(name="t", start=Coordinate("chr1", 1000),
              end=Coordinate("chr1", 1400), seq_len=300)
    run_job(job, Settings("hg38", Path("."), {}), fetch=lambda url, timeout=10: "A" * 1000)
    assert captured["PRIMER_PRODUCT_SIZE_RANGE"] == [436, 636] # 400 + 2*18, +200

def test_user_range_overrides_default(monkeypatch):
    captured = {}
    def fake(seq_args, global_args):
        captured.update(global_args)
        return {"PRIMER_PAIR_EXPLAIN": "ok"}
    monkeypatch.setattr(primer3, "design_primers", fake)

    job = Job(name="t", start=Coordinate("chr1", 1000),
              end=Coordinate("chr1", 1400), seq_len=300)
    run_job(job, Settings("hg38", Path("."), {"PRIMER_PRODUCT_SIZE_RANGE": [500, 900]}), fetch=lambda url, timeout=10: "A" * 1000)
    assert captured["PRIMER_PRODUCT_SIZE_RANGE"] == [500, 900]

@pytest.mark.parametrize(
    "option",
    [
        ["PRIMER_OPT_TM 60.5"],
        [["PRIMER_OPT_TM60.5"]],
    ]
)

def test_missing_equals(option):
     with pytest.raises(ValueError):
         assert parse_primer3_opts(option)

# Config file tests
def write_cfg(tmp_path, data):
    p = tmp_path / "cfg.json"
    p.write_text(json.dumps(data))
    return p

def test_invalid_json_rejected(tmp_path):
    p = tmp_path / "bad.json"
    p.write_text("{not_json")
    with pytest.raises(ValueError, match="not a valid JSON"):
        load_primer3_config(p)

def test_load_config_sections(tmp_path):
    p = write_cfg(tmp_path, {"global_tags": {"PRIMER_OPT_TM": 62}})
    assert load_primer3_config(p) == {"PRIMER_OPT_TM": 62}

def test_unknown_section_rejected(tmp_path):
    p = write_cfg(tmp_path, {"globals": {}})
    with pytest.raises(ValueError, match="unknown section"):
        load_primer3_config(p)

def test_global_tags_valid(tmp_path):
    p = write_cfg(tmp_path, {"global_tags": ["PRIMER_OPT_TM", 60]})
    with pytest.raises(ValueError, match="TAG: value pairs"):
        load_primer3_config(p)

def test_global_tags_not_empty(tmp_path):
    p = write_cfg(tmp_path, {"global_tags": {}})
    with pytest.raises(ValueError, match="must not be empty"):
        load_primer3_config(p)

def test_cli_overrides_file(tmp_path):
    p = write_cfg(tmp_path, {"global_tags": {"PRIMER_OPT_TM": 62, "PRIMER_MIN_SIZE": 17}})
    args = parse_args(["-s", "chr1:1000", "--p3-config", str(p),
                       "--p3-global-tags", "PRIMER_OPT_TM=58", "-o", str(tmp_path)])
    tags = settings_from_args(args).p3_global_tags
    assert tags == {"PRIMER_OPT_TM": 58, "PRIMER_MIN_SIZE": 17}

def test_missing_config(tmp_path):
    with pytest.raises(SystemExit) as exc:
        parse_args(["-s", "chr1:1000", "--p3-config", str(tmp_path / "nope.json")])
    assert exc.value.code == 2