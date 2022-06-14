from vbt3 import Molecule
from vbt3.fixed_psi import generate_dets

import sympy as sp

m = Molecule(zero_ii=True,
             subst={
                's':('S_ab','S_bc','S_cd','S_de','S_ae'),
                'h':('H_ab','H_bc','H_cd','H_de','H_ae')},
             interacting_orbs=['ab','bc','cd','de','ae']
            )

P = generate_dets(3,3,5)

res = m.Op(P[0], P[0],op='S')
