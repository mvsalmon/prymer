import unittest

target = __import__("main_new")
primer = target.Primer

class TestPrimer(unittest.TestCase):

    def test_request(self):
        test_coord = "chr1:23456789"
        parsed = primer._parse_coordinate(self, coord = test_coord)
        print(parsed)


if __name__ == '__main__':
    unittest.main()

