import pytest
from prymer import reverse_complement, build_breakpoint, breakpoint_window

## Tests for breakpoint construction

def test_reverse_complement(sample_template_4bp):
    assert reverse_complement(sample_template_4bp) == "GCAT"

@pytest.mark.parametrize(
    "rev_comp, expected",
    [
        (None, "ATGCGCAT"),
        ("start", "GCATGCAT"),
        ("end", "ATGCATGC"),
        ("both", "GCATATGC")
    ]
)
def test_build_breakpoint(sample_template_4bp, sample_template2_4bp,
                          rev_comp, expected):
    assert build_breakpoint(sample_template_4bp,
                            sample_template2_4bp,
                            rev_comp
                            ) == expected

def test_breakpoint_no_revcomp(sample_template_4bp,
                               sample_template2_4bp):
    assert build_breakpoint(sample_template_4bp,
                            sample_template2_4bp) == "ATGCGCAT"


@pytest.mark.parametrize(
    "position, expected",
    [("5", ("chr1", 500, 1000)),
     ("3", ("chr1", 1000, 1500))]
)
def test_breakpoint_window(position, expected):
    assert breakpoint_window(("chr1", 1000), position, 500) == expected