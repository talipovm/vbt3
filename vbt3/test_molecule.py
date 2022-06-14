import unittest
import sympy as sp

from vbt3 import FixedPsi, Molecule, SlaterDet
from vbt3.fixed_psi import generate_dets


class TestMolecule(unittest.TestCase):

    def test_subst_2e_1(self):
        m = Molecule(subst_2e={'J':('T_abab', 'T_bcbc')})
        self.assertEqual(
            m.subst_2e['T_abab'],
            'J'
        )

    def test_subst_2e_2(self):
        m = Molecule(subst_2e={'J':('1221', ), 'K':('1212',)})
        d1 = SlaterDet('aB')
        d2 = SlaterDet('aB')
        s = m.o2_det(d1, d2)
        self.assertEqual(
            str(sp.simplify(s)),
            '2*K'
        )

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

    def test_op_hartree_product_4(self):
        self.assertEqual(
            str(Molecule().Op_Hartree_product('AaBb', 'CcDd', op='S')),
            '(4 * S_ac*S_ac*S_bd*S_bd)'
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

    def test_op_det_fast_1(self):
        m = Molecule(zero_ii=True, subst={'s': ('S_ab', 'S_ac', 'S_bc'), 'h': ('H_ab', 'H_ac', 'H_bc')})
        m.generate_basis(2, 2, 3)
        s1 = m.basis[0].determinants[0]['det_string']
        sd = SlaterDet(s1)
        m = m.op_det(sd, sd)
        self.assertEqual(
            str(m),
            '-2*h*s*(2 - 2*s**2)'
        )

    def test_op_det_fast_2(self):
        m = Molecule(zero_ii=True, subst={'s': ('S_ab', 'S_ac', 'S_bc'), 'h': ('H_ab', 'H_ac', 'H_bc')})
        m.generate_basis(2, 2, 3)
        s1 = m.basis[0].determinants[0]['det_string']
        sd = SlaterDet(s1)
        m = m.op_det(sd, sd, op='S')
        self.assertEqual(
            str(m),
            '(2 - 2*s**2)**2'
        )

    def test_op_det_fast_3(self):
        m = Molecule(zero_ii=True,
                     subst={
                         's': ('S_ab', 'S_bc', 'S_cd'),
                         'h': ('H_ab', 'H_bc', 'H_cd')},
                     interacting_orbs=['ab', 'bc', 'cd']
                     )
        m.generate_basis(2, 2, 4)
        M2 = m.build_matrix(m.basis, op='S')
        self.assertEqual(
            str(M2[0,1]),
            '2*s*(2 - 2*s**2)'
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
            'Matrix([[4*s**2 + 4, 4*s**2 + 4*s, 4*s**2 + 4*s], [4*s**2 + 4*s, 4*s**2 + 4, 4*s**2 + 4], [4*s**2 + 4*s, 4*s**2 + 4, 4*s**2 + 4]])'
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

    def test_o2_1(self):
        m = Molecule(zero_ii=True, subst={'s': ('S_ab',), 'h': ('H_ab',)}, interacting_orbs=['ab', ])
        d1 = SlaterDet('aB')
        d2 = SlaterDet('aB')
        s = m.o2_det(d1, d2)
        self.assertEqual(
            str(sp.simplify(s)),
            '2*T_abab'
        )

    def test_o2_2(self):
        m = Molecule(zero_ii=True, subst={'s': ('S_ab',), 'h': ('H_ab',)}, interacting_orbs=['ab', ])
        d1 = SlaterDet('ab')
        d2 = SlaterDet('ab')
        s = m.o2_det(d1, d2)
        self.assertEqual(
            str(sp.simplify(s)),
            '-2*T_aabb + 2*T_abab'
        )

    def test_o2_3(self):
        m = Molecule(zero_ii=True,
                     subst={'s': ('S_ab', 'S_bc', 'S_cd'),
                            'h': ('H_ab', 'H_bc', "H_cd")},
                     interacting_orbs=['ab', 'bc', 'cd'])
        d1 = SlaterDet('aBc')
        d2 = SlaterDet('aBd')
        s = m.o2_det(d1, d2)
        self.assertEqual(
            str(sp.simplify(s)),
            '-6*T_aadc + 6*T_abab*s + 6*T_acad + 6*T_bcbd'
        )


if __name__ == '__main__':
    unittest.main()


