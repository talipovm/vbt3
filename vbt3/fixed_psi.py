import vbt3
from vbt3.functions import generate_det_strings
from vbt3.functions import attempt_int


class FixedPsi:
    # Perfect-pair expansion of determinants
    def __init__(self, x=None, coupled_pairs=None):
        """

        """
        self.Nel = 0

        self.dets = []
        self.coefs = []

        if x is None:
            return
        elif x.__class__.__name__ == 'str':
            self.add_str_det(x)
        elif x.__class__.__name__ == 'SlaterDet':
            self.add_det(x)
        elif x.__class__.__name__ == 'FixedPsi':
            self.Nel = x.Nel
            self.dets = x.dets
            self.coefs = x.coefs

        if coupled_pairs is not None:
            for i, j in coupled_pairs:
                self.couple_orbitals(i, j)

    def __getitem__(self, item):
        return self.dets[item]

    def __len__(self):
        return len(self.dets)

    def __contains__(self, item):
        for d in self.dets:
            if d.det_string == item:
                return True
        return False

    def contains_det(self, det_string):
        for D in self.dets:
            if det_string == D['det_string']:
                return True
        return False

    def add_det(self, det, coef=+1):
        assert det.__class__.__name__ == 'SlaterDet'
        assert det not in self
        assert self.Nel in (0, det.Nel)

        self.dets.append(det)
        cf = attempt_int(coef)
        self.coefs.append(cf)
        if self.Nel == 0:
            self.Nel = det.Nel

    def add_str_det(self, det_string, coef=+1):
        sd = vbt3.SlaterDet(det_string)
        self.add_det(sd)

    def add_fixedpsi(self, p, coef=1.0):
        for i in range(len(p)):
            self.add_det(p.dets[i], p.coefs[i] * coef)

    def couple_orbitals(self, o1, o2):
        # generate determinants that represent a singlet bonding coupling between two orbitals.
        # Orbital numbering starts from 0
        determinants = self.dets.copy()
        coefs = self.coefs.copy()
        for i in range(len(determinants)):
            ds = determinants[i].det_string
            coef = coefs[i]
            # Flip spins
            c1, c2 = [c.lower() if c.isupper() else c.upper() for c in [ds[o1], ds[o2]]]
            assert c1.lower() != c2.lower(), 'Cannot couple the same orbital'
            # Flip positions
            ds2 = ds[:o1] + c2 + ds[(o1 + 1):o2] + c1 + ds[(o2 + 1):]
            self.add_str_det(ds2, coef=coef)

    def __repr__(self):
        s = ''
        for i in range(len(self)):
            d = self.dets[i]
            cf = self.coefs[i]
            dc = attempt_int(cf)

            if dc > 0:
                if dc == 1.0:
                    s += '+'
                else:
                    s += '+%s' % dc
            elif dc < 0:
                if dc == -1.0:
                    s += '-'
                else:
                    s += '-%s' % dc
            s += str(d)
        if s[0] == '+':
            s = s[1:]
        return s


def generate_dets(Nela, Nelb, Norb):
    """
    Generate all possible determinants for a given number of electrons and atomic orbitals.
    :param Nela: Number of alpha electrons
    :param Nelb: Number of beta electrons
    :param Norb: Number of atomic orbitals
    :return: List of FixedPsi objects, each containing one determinant
    """
    L = generate_det_strings(Nela, Nelb, Norb)
    PP = [None,]*len(L)
    for i in range(len(L)):
        PP[i] = FixedPsi(L[i])
    return PP
