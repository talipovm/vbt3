from vbt3.fixed_psi import generate_dets
from vbt3 import Molecule


m = Molecule(zero_ii=True,
             subst={
                's':('S_ab', 'S_bc', 'S_cd'),
                'h':('H_ab', 'H_bc', 'H_cd')},
             interacting_orbs=['ab', 'bc', 'cd']
            )

P = generate_dets(2,2,4)

M = m.build_matrix(P, op='S')
