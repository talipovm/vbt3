import unittest
from vbt3 import SlaterDet


class TestSlaterDet(unittest.TestCase):

    def test_add_1(self):
        self.assertEqual(
            str(SlaterDet('aB') + 'bA'),
            '|aB|+|bA|'
        )

    def test_SlaterDet_1(self):
        self.assertEqual(
            SlaterDet('AbCd').spins,
            '-+-+'
        )

    def test_SlaterDet_2(self):
        self.assertEqual(
            SlaterDet('AbCd').det_string,
            'AbCd'
        )

    def test_SlaterDet_3(self):
        self.assertEqual(
            SlaterDet('AbCd').alpha_indices,
            [1, 3]
        )

    def test_SlaterDet_4(self):
        self.assertEqual(
            SlaterDet('AbCd').get_orbital_permutations(),
            [['AbCd', 'CbAd', 'AdCb', 'CdAb'], [1, -1, -1, 1]]
        )

    def test_SlaterDet_5(self):

        self.assertEqual(
            SlaterDet('').get_orbital_permutations(),
            [[''], [1]]
        )

    def test_SlaterDet_6(self):
        self.assertEqual(
            SlaterDet('A').get_orbital_permutations(),
            [['A'], [1]]
        )

    def test_SlaterDet_7(self):
        self.assertEqual(
            str(SlaterDet('aBcD')),
            '|aBcD|'
        )

if __name__ == '__main__':
    unittest.main()