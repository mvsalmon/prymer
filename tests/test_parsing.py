import pytest
from prymer import parse_coordinates

def test_parse_simple_coordinates():
    assert parse_coordinates("chr1:123456") == ("chr1", 123456)

def test_parse_strip_commas():
    assert parse_coordinates("chr1:12,3456") == ("chr1", 123456)

def test_parse_rejects_bad_format():
    with pytest.raises(ValueError, match="Invalid coordinates"):
        parse_coordinates("1:123456") # missing 'chr'

@pytest.mark.parametrize(
    "text, expected",
    [
        ("chr1:1000", ("chr1", 1000)),
        ("chr22:5", ("chr22", 5)),
        ("chrX:7000", ("chrX", 7000)),
        ("chrY:7000", ("chrY", 7000)),
        ("chr10:1,000,000", ("chr10", 1_000_000)),
    ],
)
def test_parse_valid(text, expected):
    assert parse_coordinates(text) == expected

@pytest.mark.parametrize(
    "bad",
    ["1:100", "chr23:10000", "chr0:100000", "chr1:", "chr:100000", "chr1-1000", "chr1:1.000.000", "chrM:1000",
     " chr1:1000"],
)
def test_parse_invalid(bad):
    with pytest.raises(ValueError):
        parse_coordinates(bad)