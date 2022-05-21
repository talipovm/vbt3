import unittest

from vbt3.fixed_psi import attempt_int


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