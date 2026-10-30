import pytest
from prymer import parse_primer3_opts

def test_comma_range():
    result = parse_primer3_opts(["PRIMER_PRODUCT_SIZE_RANGE=100-300,300-500"])
    assert isinstance(result["PRIMER_PRODUCT_SIZE_RANGE"], list)

@pytest.mark.parametrize(
    "options, expected",
    [
        (["PRIMER_OPT_TM=60.5"],  {"PRIMER_OPT_TM": 60.5}),
        (["PRIMER_MIN_THREE_PRIME_DISTANCE=-1"],  {"PRIMER_MIN_THREE_PRIME_DISTANCE": -1}),
        (["SEQUENCE_INCLUDED_REGION=50,400"],  {"SEQUENCE_INCLUDED_REGION": [50, 400]}),
        (["SEQUENCE_TARGET=50,2 100,2"],   {"SEQUENCE_TARGET": [[50, 2], [100, 2]]}),
        (["PRIMER_PRODUCT_SIZE_RANGE=100,300 400,500"],
         {"PRIMER_PRODUCT_SIZE_RANGE": [[100, 300], [400, 500]]}),
        (["SEQUENCE_ID=sample.1"], {"SEQUENCE_ID": 'sample.1'})
    ],
)

def test_p3_option_parsing(options, expected):
    assert parse_primer3_opts(options) == expected