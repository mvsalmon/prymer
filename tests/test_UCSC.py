import pytest
from prymer import build_UCSC_url, reverse_complement, build_breakpoint
from tests.conftest import sample_template2_4bp


def test_url_construction():
    assert (build_UCSC_url("hg38", "chr1", 10, 20) ==
            "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr1;start=10;end=20")

def test_reverse_complement(sample_template_4bp):
    assert reverse_complement(sample_template_4bp) == "GCAT"

