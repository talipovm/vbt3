import unittest
from vbt3 import SlaterDet


class TestSlaterDet(unittest.TestCase):
    def test_SlaterDet_1(self):
        self.assertEqual(
            SlaterDet('AbCd').spins,
            '+-+-'
        )

    def test_SlaterDet_2(self):
        self.assertEqual(
            SlaterDet('AbCd').det_string,
            'AbCd'
        )

    def test_SlaterDet_3(self):
        self.assertEqual(
            SlaterDet('AbCd').alpha_indices,
            [0, 2]
        )

    def test_SlaterDet_4(self):
        self.assertEqual(
            SlaterDet('AbCd').get_orbital_permutations(),
            [['AbCd', 'AdCb', 'CbAd', 'CdAb'], [1, -1, -1, 1]]
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

if __name__ == '__main__':
    unittest.main()