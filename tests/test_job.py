import pytest
from pathlib import Path
from prymer import run_job, Job, Coordinate, Settings, build_UCSC_url


def make_fake_fetch(responses: dict[str, str]):
    def fake(url, timeout=10):
        return responses[url]
    return fake

def test_run_job_fusion(dummy_response_sequence):
    seq = dummy_response_sequence.upper()
    job = Job(name="test", start=Coordinate("chr1", 1000),
              end=Coordinate("chr1", 5000), sv= True, seq_len=500,
              end_primer_position="3")
    settings = Settings("hg38", Path("."), {})
    fake = make_fake_fetch({
        build_UCSC_url("hg38", "chr1", 500, 1000): seq[2500:3000],
        build_UCSC_url("hg38", "chr1", 5000, 5500): seq[3500:4000],
    })
    result = run_job(job, settings, fetch=fake)
    assert result.job is job
    assert result.explain