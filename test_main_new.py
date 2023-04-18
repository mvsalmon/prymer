from unittest import TestCase

from main_new import Primer

class TestPrimer(TestCase):

    def test_ucsc_request_single(self):
        """Test if UCSC API request returns a dict when using a single coordinate"""
        p = Primer(start_coordinate="chr1:123456789")
        assert isinstance(p.UCSC_response, dict)
    def test_ucsc_request_double(self):
        """Test if UCSC API request returns a dict when using two coordinates (non-fusion)"""
        p = Primer(start_coordinate="chr1:123456789",
                   end_coordinate="chr1:123456899")
        assert isinstance(p.UCSC_response, dict)

    def test_ucsc_request_breakpoint_start(self):
        """Test if UCSC API request returns a dict when 2 coordinates on different chroms and
        --fusion_breakpoint specified"""
        p = Primer(start_coordinate="chr1:123456789",
                   end_coordinate="chr5:123456899",
                   fusion_breakpoint=True)
        assert isinstance(p.UCSC_start_breakpoint_response, dict)

    def test_ucsc_request_breakpoint_end(self):
        """Test if UCSC API request returns a dict when 2 coordinates on different chroms and
        --fusion_breakpoint specified"""
        p = Primer(start_coordinate="chr1:123456789",
                   end_coordinate="chr5:123456789",
                   fusion_breakpoint=True)
        assert isinstance(p.UCSC_end_breakpoint_response, dict)

