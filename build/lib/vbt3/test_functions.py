from unittest import TestCase
from vbt3.functions import generate_det_strings, attempt_int


class TestFixedPsi(TestCase):

    def test_generate_det_strings(self):
        self.assertEqual(
            str(generate_det_strings(2, 2, 3)),
            "['aAbB', 'aAbC', 'aBbC', 'aAcB', 'aAcC', 'aBcC', 'bAcB', 'bAcC', 'bBcC']"
        )

