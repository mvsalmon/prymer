import pytest

from prymer.cli import parse_args
from prymer.snps import load_maf_manifest

def test_read_sources(tmp_path):
    path = tmp_path / "sources.txt"
    path.write_text(
        "# Order of sources for MAF values found in dbSNP155\n"
        "1000Genomes\n"
        "dbGaP_PopFreq\n"
        "   # indented comment line\n"
        "TOPMED\n"
        "KOREAN\n"
        "SGDP_PRJ\n",
   )
    assert load_maf_manifest(path) == [
        "1000Genomes",
        "dbGaP_PopFreq",
        "TOPMED",
        "KOREAN",
        "SGDP_PRJ"
    ]

def test_maf_manifest_matches_source_count():
    assert len(load_maf_manifest()) == 31

def test_maf_manifest_unique():
    sources = load_maf_manifest()
    assert len(sources) == len(set(sources))

def test_maf_manifest_gnomad_position():
    assert load_maf_manifest().index("GnomAD") == 12

def test_manifest_not_loaded_without_snp_track():
    args = parse_args(["-s", "chr1:10000"])
    assert args.maf_manifest is None

def test_manifest_not_loaded_for_any_source():
    args = parse_args(["-s", "chr1:10000", "--snp-track", "dbSnp155"])
    assert args.maf_manifest is None
    assert args.maf_source is None
    assert args.snp_freq_cutoff == 0.01

def test_manifest_loaded_for_named_source():
    args = parse_args(["-s", "chr1:10000", "--snp-track", "dbSnp155",
                       "--maf-source", "GnomAD", "TOPMED"])
    assert "GnomAD" in args.maf_manifest

def test_cutoff_default_applies_with_maf_source():
    args = parse_args(["-s", "chr1:10000", "--snp-track", "dbSnp155",
                       "--maf-source", "GnomAD"])
    assert args.snp_freq_cutoff == 0.01