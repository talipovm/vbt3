import string
from vbt3 import FixedPsi
from itertools import combinations


def generate_det_strings(Na, Nb, Norbs):
    """
    Generate all possible determinant strings for a given number of electrons and atomic orbitals.
    :param Na: Number of alpha electrons
    :param Nb: Number of beta electrons
    :param Norbs: Number of atomic orbitals
    :return: list of determinant strings
    """
    assert Na >= Nb, 'Na cannot be smaller than Nb'
    result = []
    for a in combinations(string.ascii_lowercase[:Norbs], Na):
        for b in combinations(string.ascii_uppercase[:Norbs], Nb):
            s = ''
            for i in range(Nb):
                s += a[i] + b[i]
            for i in range(Nb,Na):
                s += a[i]
            result.append(s)
    return result


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