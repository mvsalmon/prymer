import pytest
from pathlib import Path
from prymer import (run_job, Job, Coordinate, Settings, Result, parse_args, jobs_from_args, build_UCSC_url,
                    write_output, main)


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

def test_parse_args_builds_coorindates():
    args = parse_args(["-s", "chr1:10000"])
    assert args.start_coordinate == Coordinate("chr1", 10000)

def test_fusion_requires_end(capsys):
    with pytest.raises(SystemExit) as exc:
        parse_args(["-s", "chr1:10000", "-f"])
        assert exc.value.code == 2
        assert "--end-coordinate" in capsys.readouterr().err

def test_default_name_uses_both_coordinates():
        args = parse_args(["-s", "chr1:10000", "-e", "chr6:54321", "-f"])
        assert jobs_from_args(args)[0].name == "chr1_10000_chr6_54321"

def test_different_chromosomes_without_fusion_rejected():
    with pytest.raises(ValueError, match="differ"):
        Job(name="x", start=Coordinate("chr1", 10000), end=Coordinate("chr2", 20000))

def test_write_output_names_file_after_job(tmp_path):
    job = Job(name="demo", start=Coordinate("chr1", 1000))
    result = Result(job, pairs={"PAIR_0": {"PRIMER_PAIR_PENALTY": 0.1}}, explain="ok")
    assert write_output(result, tmp_path) == tmp_path / "demo.csv"

def test_write_output_skips_empty(tmp_path):
    job = Job(name="demo", start=Coordinate("chr1", 1000))
    assert write_output(Result(job, {}, "considered 0"), tmp_path) is None
    assert list(tmp_path.iterdir()) == []

def test_main_returns_1_for_missing_output_dir(tmp_path):
    missing = tmp_path / "nope"
    assert main(["-s", "chr1:1000", "-o", str(missing)]) == 1

def test_negative_start():
    with pytest.raises(ValueError):
        Job(name="x", start=Coordinate("chr1", 200))

def test_end_before_start():
    with pytest.raises(ValueError):
        Job(name="x", start=Coordinate("chr1", 10000), end=Coordinate("chr1", 9500))