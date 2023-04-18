from unittest import TestCase

from main_new import Primer

class TestPrimer(TestCase):

    def test_ucsc_request_single(self):
        p = Primer(start_coordinate="chr1:123456789")
        assert isinstance(p.UCSC_response, dict)

    def test_ucsc_request_double(self):
        p = Primer(start_coordinate="chr1:123456789",
                   end_coordinate="chr1:123456899")
        assert isinstance(p.UCSC_response, dict)
