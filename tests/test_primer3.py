import pytest
from prymer import parse_primer3_opts


def test_comma_range():
    result = parse_primer3_opts(["PRIMER_PRODUCT_SIZE_RANGE=100-300,300-500"])
    assert isinstance(result["PRIMER_PRODUCT_SIZE_RANGE"], list)