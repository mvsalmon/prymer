import pytest
from prymer.ucsc import build_ucsc_url, fetch_sequence


## Tests for UCSC related functions

class FakeResponse:
    def __init__(self, payload): self._payload = payload
    def raise_for_status(self): pass
    def json(self): return self._payload

def test_fetch_sequence_uppercases(monkeypatch):
    monkeypatch.setattr("requests.get",
                        lambda url, timeout: FakeResponse({"dna": "atgcNNac"}))
    assert fetch_sequence("http://x") == "ATGCNNAC"

@pytest.mark.network
def test_ucsc_live(dummy_response_sequence):
    url = build_ucsc_url("hg38", "chr1",
                         10000,
                         15000)
    assert fetch_sequence(url) == dummy_response_sequence.upper()

def test_url_construction():
    assert (build_ucsc_url("hg38", "chr1", 10, 20) ==
            "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr1;start=10;end=20")