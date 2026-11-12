import unittest

from hplc_autoquant.single_point_assay import calculate_single_point


class TestSinglePointAssay(unittest.TestCase):
    def test_external_standard_ratio(self):
        result = calculate_single_point(1000, 0.5, 250, 2)
        self.assertAlmostEqual(result['unknown_concentration'], 0.25)

    def test_response_factor(self):
        result = calculate_single_point(1000, 1.0, 500, response_factor=2)
        self.assertAlmostEqual(result['unknown_concentration'], 0.25)

    def test_rejects_invalid_area(self):
        with self.assertRaises(ValueError):
            calculate_single_point(0, 1, 10)
