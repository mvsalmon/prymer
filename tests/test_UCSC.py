import pytest
from prymer import build_UCSC_url, fetch_sequence


## Tests for UCSC related functions

def test_url_construction():
    assert (build_UCSC_url("hg38", "chr1", 10, 20) ==
            "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr1;start=10;end=20")


def test_ucsc_call(dummy_response_sequence):
    url = build_UCSC_url("hg38", "chr1",
                         10000,
                         15000)
    seq = fetch_sequence(url)
    assert seq == dummy_response_sequence
