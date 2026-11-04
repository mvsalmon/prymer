import pytest
import copy
import re
from prymer import parse_primer3_opts, parse_primer3_output

def test_comma_range():
    result = parse_primer3_opts(["PRIMER_PRODUCT_SIZE_RANGE=100-300,300-500"])
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

def test_parse_primer3_output(p3_raw_output):
    result = parse_primer3_output(p3_raw_output)
    assert "PAIR_0" in result