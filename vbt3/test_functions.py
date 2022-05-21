from unittest import TestCase
from vbt3.functions import generate_det_strings, attempt_int


class TestFixedPsi(TestCase):

    def test_generate_det_strings(self):
        self.assertEqual(
            str(generate_det_strings(2, 2, 3)),
            "['aAbB', 'aAbC', 'aBbC', 'aAcB', 'aAcC', 'aBcC', 'bAcB', 'bAcC', 'bBcC']"
        )

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
