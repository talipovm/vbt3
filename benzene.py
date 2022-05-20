import pickle
from vbt3 import Molecule

m = Molecule(zero_ii=True,
             subst={
                's':('S_ab','S_bc','S_cd','S_de','S_ef','S_af'),
                'h':('H_ab','H_bc','H_cd','H_de','H_ef','H_af')},
             interacting_orbs=['ab','bc','cd','de','ef','af']
            )

dir = 'C:\\Users\\talip\\OneDrive - New Mexico State University\\NMSU-LAPTOP-SGUG66BT\Research\Jupiter\\Valence Bond Theory\\'
fname = dir + 'benzene.pickle'

file = open(fname, 'rb')
benzene = pickle.load(file)
file.close()

print('Started calculating the energy')
E = m.energy(benzene)
print('Done, saving the results')

fname = dir + 'benzene-energy.pickle'
file = open(fname, 'wb')
pickle.dump(E, file)
file.close()