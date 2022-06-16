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
            fp[0].det_string,
            'abc'
        )

    def test_get_item_2(self):
        fp = FixedPsi('abc')
        fp.add_str_det('bcd')
        self.assertEqual(
            fp[1].det_string,
            'bcd'
        )

    def test_len(self):
        fp = FixedPsi('abc')
        fp.add_str_det('bcd')
        self.assertEqual(
            len(fp),
            2
        )

    def test_add_str_1(self):
        fp2 = FixedPsi('aA')
        fp2.add_str_det('aB', coef=-1)
        fp2.add_str_det('bA', coef=-1)
        fp2.add_str_det('bB')
        self.assertEqual(
            str(fp2),
            '|aA|-|aB|-|bA|+|bB|'
        )

    def test_iter_1(self):
        fp = FixedPsi('aA')
        fp.add_str_det('bB')
        s = ''
        for d,c in fp:
            s += '%s: %s;  ' % (d,c)
        self.assertEqual(
            s,
            '|aA|: 1;  |bB|: 1;  '
        )
