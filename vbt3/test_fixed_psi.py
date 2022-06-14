import unittest

from vbt3.fixed_psi import attempt_int, FixedPsi


class TestMolecule(unittest.TestCase):
    def test_attempt_int_1(self):
        self.assertEqual(
            str(attempt_int(3.5)),
            '3.5'
        )

    def test_attempt_int(self):
        self.assertEqual(
            str(attempt_int(3.0)),
            '3'
        )

    def test_get_item_1(self):
        fp = FixedPsi('abc')
        self.assertEqual(
            fp[0]['det_string'],
            'abc'
        )

    def test_get_item_2(self):
        fp = FixedPsi('abc')
        fp.add_str_det('bcd')
        self.assertEqual(
            fp[1]['det_string'],
            'bcd'
        )

    def test_len(self):
        fp = FixedPsi('abc')
        fp.add_str_det('bcd')
        self.assertEqual(
            len(fp),
            2
        )
