import unittest

from vbt3 import FixedPsi, Molecule
from vbt3.fixed_psi import generate_dets


class TestMolecule(unittest.TestCase):

    def test_op_hartree_product_1(self):
        self.assertEqual(
            str(Molecule().Op_Hartree_product('AaBb', 'CcDd')),
            '(H_ac*S_ac*S_bd*S_bd + S_ac*H_ac*S_bd*S_bd + S_ac*S_ac*H_bd*S_bd + S_ac*S_ac*S_bd*H_bd)'
        )

    def test_op_hartree_product_2(self):
        self.assertEqual(
            str(Molecule().Op_Hartree_product('aAbBcC', 'aAbBcC', op='H')),
            '(0)'
        )

    def test_op_hartree_product_3(self):
        self.assertEqual(
            str(Molecule().Op_Hartree_product('aA', 'aB', op='H')),
            '(H_ab)'
        )
    def test_op_1(self):
        self.assertEqual(
            str(Molecule(zero_ii=False).Ops('AbCd', 'AbCd')),
            '-H_aa*S_bd**2 + H_aa + 2*H_ac*S_ac*S_bd**2 - 2*H_ac*S_ac - H_bb*S_ac**2 + H_bb + 2*H_bd*S_ac**2*S_bd - 2*H_bd*S_bd - H_cc*S_bd**2 + H_cc - H_dd*S_ac**2 + H_dd'
        )

    def test_op_2(self):
        self.assertEqual(
            str(Molecule().Ops('', '')),
            '0'
        )

    def test_op_3(self):
        self.assertEqual(
            str(Molecule().Ops('', '', op='S')),
            '1'
        )

    def test_op_4(self):
        self.assertEqual(
            str(Molecule(zero_ii=False).Ops('A', 'A')),
            'H_aa'
        )

    def test_op_5(self):
        self.assertEqual(
            str(Molecule(zero_ii=False).Ops('b', 'b')),
            'H_bb'
        )

    def test_op_6(self):
        self.assertEqual(
            str(Molecule(zero_ii=False).Ops('AA', 'BB')),
            '0'
        )

    def test_op_7(self):
        self.assertEqual(
            str(Molecule(zero_ii=False, interacting_orbs=['ab']).Ops('AaBb', 'CcDd')),
            '0'
        )

    def test_op_8(self):
        self.assertEqual(
            str(Molecule(zero_ii=False, interacting_orbs=['ab', 'bc', 'cd']).Ops('ABbCcd', 'ABbCcd', op='S')),
            '6*(S_ab**2 + S_bc**2 - 1)*(S_bc**2 + S_cd**2 - 1)'
        )

    def test_build_matrix(self):
        m = Molecule(zero_ii=True, subst={'s': ('S_ab', 'S_ac', 'S_bc'), 'h': ('H_ab', 'H_ac', 'H_bc')})
        P = [0, ] * 3
        P[0] = FixedPsi('Ab', coupled_pairs=[(0, 1)])
        P[1] = FixedPsi('Bc', coupled_pairs=[(0, 1)])
        P[2] = FixedPsi('Bc', coupled_pairs=[(0, 1)])
        m = m.build_matrix(P, op='S')
        self.assertEqual(
            str(m),
            'Matrix([[4*(s**2 + 1), 4*s*(s + 1), 4*s*(s + 1)], [4*s*(s + 1), 4*(s**2 + 1), 4*(s**2 + 1)], [4*s*(s + 1), 4*(s**2 + 1), 4*(s**2 + 1)]])'
        )

    def test_normalized_self_energy(self):
        m = Molecule(zero_ii=False)
        P = FixedPsi('AB')
        self.assertEqual(
            str(m.energy(P)),
            '-(2*H_aa - 4*H_ab*S_ab + 2*H_bb)/(2*(S_ab - 1)*(S_ab + 1))'
        )

    def test_couple(self):
        m = Molecule(zero_ii=True, subst={'s': ('S_ab', 'S_ac', 'S_bc'), 'h': ('H_ab', 'H_ac', 'H_bc')})
        P = generate_dets(1, 1, 3)
        result = m.couple(P)
        self.assertEqual(
            str(result),
            '[|aA|+|aB|+|aC|+|bA|+|bB|+|bC|+|cA|+|cB|+|cC|]'
        )

    def test_couple2(self):
        m = Molecule(zero_ii=True, subst={'s': ('S_ab', 'S_bc'), 'h': ('H_ab', 'H_bc')}, interacting_orbs=['ab', 'bc'])
        P = generate_dets(1, 1, 3)
        result = m.couple(P, N_tries=20)
        self.assertEqual(
            str(result),
            '[|aA|+|aC|+2|bB|+|cA|+|cC|, |aB|+|bA|+|bC|+|cB|]'
        )


if __name__ == '__main__':
    unittest.main()

