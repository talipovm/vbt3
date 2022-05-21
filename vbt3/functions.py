import string
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


def attempt_int(x):
    """
    Convert a float to int if it is numerically equivalent
    Parameters
    ----------
    x: a number (float)
    Returns
    -------
    int or float
    """
    if x == int(x):
        return int(x)
    else:
        return x
