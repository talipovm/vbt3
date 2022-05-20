

class FixedPsi:
    # Perfect-pair expansion of determinants
    def __init__(self, X=None, coupled_pairs=[]):
        """

        """
        self.determinants = []
        self.Nel = 0

        if X is None:
            return

        if X.__class__.__name__ == 'str':
            self.add_str_det(X)

        if X.__class__.__name__ == 'SlaterDet':
            self.add_str_det(X.det_string)

        if X.__class__.__name__ == 'FixedPsi':
            self.determinants = X.determinants
            self.Nel = X.Nel

        for i, j in coupled_pairs:
            self.couple_orbitals(i, j)

    def contains_det(self, det_string):
        for D in self.determinants:
            if det_string == D['det_string']:
                return True
        return False

    def add_str_det(self, det_string, coef=+1):
        if not isinstance(det_string, str):
            raise Exception('ds is %s instead of str' % (type(det_string)))
        if self.contains_det(det_string):
            raise Exception('New determinant %s is already in the list' % (det_string))
        self.determinants.append({'det_string':det_string,'coef':coef})
        self.Nel = len(det_string)

    def add_FixedPsi(self, p, coef=1.0):
        for d in p.determinants:
            self.add_str_det(d['det_string'], d['coef'] * coef)

    def couple_orbitals(self, o1, o2):
        # generate determinants that represent a singlet bonding coupling between two orbitals.
        # Orbital numbering starts from 0
        determinants = self.determinants.copy()
        for D in determinants:
            # Flip spins
            ds = D['det_string']
            c1, c2 = [c.lower() if c.isupper() else c.upper() for c in [ds[o1], ds[o2]]]
            if c1.lower() == c2.lower():
                raise Exception('Cannot couple the same orbital, (%s,%s)' % (c1, c2))
            # Flip positions
            ds2 = ds[:o1] + c2 + ds[(o1 + 1):o2] + c1 + ds[(o2 + 1):]
            self.add_str_det(ds2, coef=+D['coef'])

    def __repr__(self):
        s = ''
        for d in self.determinants:
            dc = d['coef']
            if int(dc)==dc:
                dc = int(dc)

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
            s += '|%s|' % d['det_string']
        if s[0] == '+':
            s = s[1:]
        return s

