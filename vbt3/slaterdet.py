from vbt3.orbital_permutations import OrbitalPermutations
import vbt3
import logging

logging.basicConfig(format=('%(levelname)-8s: %(message)s'))


class SlaterDet:
    """
    Parses the Slater determinant string |s| and computes the <det1|O|det2> elements
    """
    def __init__(self, s=""):
        self.det_string = s
        self.alpha_indices = []
        self.alpha_string = ''
        self.beta_indices = []
        self.beta_string = ''
        self.spins = ''
        self.Nel = 0
        if len(self.det_string) > 0:
            self.parse_det()

    def __repr__(self):
        s = '|%s|' % self.det_string
        return s

    def __add__(self, other):
        result = vbt3.FixedPsi(self)
        result.add_str_det(other)
        return result

    def parse_det(self):
        s = self.det_string
        if len(s) == 0:
            return

        i = 0
        for c in s:
            if c.islower():
                self.alpha_indices.append(i)
                self.spins += '+'
                self.alpha_string += c
            else:
                self.beta_indices.append(i)
                self.spins += '-'
                self.beta_string += c
            i = i + 1
        self.Nel = i

    def get_orbital_permutations(self):
        # gets all spin-restricted permutations of orbital products
        A = OrbitalPermutations(len(self.alpha_indices))
        B = OrbitalPermutations(len(self.beta_indices))

        dets = []
        signs = []
        for a_orbs, a_sign in zip(A.permutations, A.permutation_signs):
            for b_orbs, b_sign in zip(B.permutations, B.permutation_signs):
                i_a = 0
                i_b = 0
                s = ''
                for i in range(len(self.det_string)):
                    if self.spins[i] == '+':
                        c = self.alpha_string[a_orbs[i_a]]
                        i_a = i_a + 1
                    else:
                        c = self.beta_string[b_orbs[i_b]]
                        i_b = i_b + 1
                    s = s + c
                dets.append(s)
                signs.append(a_sign * b_sign)

        return ([dets, signs])

    def is_compatible(self, R):
        if self.Nel != R.Nel:
            logging.warning('Different number of electrons: %i vs %i' % (self.Nel, R.Nel))
            return (False)
        if self.spins != R.spins:
            logging.warning('The determinant spins are incompatible: %s vs %s' % (self.spins, R.spins))
            return (False)
        return (True)


if __name__ == '__main__':
    print(SlaterDet('AbCd').spins)
