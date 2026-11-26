import json, pytest
from prymer.p3 import load_primer3_config
from prymer.cli import parse_args, settings_from_args

def write_cfg(tmp_path, data):
    p = tmp_path / "cfg.json"
    p.write_text(json.dumps(data))
    return p

def test_invalid_json_rejected(tmp_path):
    p = tmp_path / "bad.json"
    p.write_text("{not_json")
    with pytest.raises(ValueError, match="not a valid JSON"):
        load_primer3_config(p)

def test_load_config_sections(tmp_path):
    p = write_cfg(tmp_path, {"global_tags": {"PRIMER_OPT_TM": 62}})
    assert load_primer3_config(p) == {"PRIMER_OPT_TM": 62}

def test_unknown_section_rejected(tmp_path):
    p = write_cfg(tmp_path, {"globals": {}})
    with pytest.raises(ValueError, match="unknown section"):
        load_primer3_config(p)

def test_global_tags_valid(tmp_path):
    p = write_cfg(tmp_path, {"global_tags": ["PRIMER_OPT_TM", 60]})
    with pytest.raises(ValueError, match="TAG: value pairs"):
        load_primer3_config(p)

def test_global_tags_not_empty(tmp_path):
    p = write_cfg(tmp_path, {"global_tags": {}})
    with pytest.raises(ValueError, match="must not be empty"):
        load_primer3_config(p)

def test_cli_overrides_file(tmp_path):
    p = write_cfg(tmp_path, {"global_tags": {"PRIMER_OPT_TM": 62, "PRIMER_MIN_SIZE": 17}})
    args = parse_args(["-s", "chr1:1000", "--p3-config", str(p),
                       "--p3-global-tags", "PRIMER_OPT_TM=58", "-o", str(tmp_path)])
    tags = settings_from_args(args).p3_global_tags
    assert tags == {"PRIMER_OPT_TM": 58, "PRIMER_MIN_SIZE": 17}

def test_missing_config(tmp_path):
    with pytest.raises(SystemExit) as exc:
        parse_args(["-s", "chr1:1000", "--p3-config", str(tmp_path / "nope.json")])
    assert exc.value.code == 2