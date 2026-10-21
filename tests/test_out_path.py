import os
import sys
import pytest
from prymer import validate_output_path

def test_validiate_output_path_ok(tmp_path):
    assert validate_output_path(tmp_path) == tmp_path

def test_validiate_output_path_missing(tmp_path):
    with pytest.raises(FileNotFoundError):
        validate_output_path(str(tmp_path / "does_not_exist"))

# Test if writable
@pytest.mark.skipif(sys.platform == "win32" or os.geteuid() == 0,
                    reason="permision bits unreliable on Windows/as root")
def test_validate_output_path_not_writable(tmp_path):
    tmp_path.chmod(0o500) # r+x only
    try:
        with pytest.raises(PermissionError):
            validate_output_path(str(tmp_path))
    finally:
        tmp_path.chmod(0o700) # restore permission for cleanup