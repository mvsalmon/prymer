import pytest
import requests
from prymer import build_UCSC_url, reverse_complement, build_breakpoint

## Tests for USCS related functions

def test_url_construction():
    assert (build_UCSC_url("hg38", "chr1", 10, 20) ==
            "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr1;start=10;end=20")


def test_ucsc_call()
