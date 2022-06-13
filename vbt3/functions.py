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


def place_low(s, i):
    """
    Flips two characters in the string to place the next available lower-case character at i-th position.
    E.g. place_low('aBcD',1) will return 'acBD'
    Parameters
    ----------
    s: a string containing upper and lower case characters
    i: position in the string

    Returns
    -------
    a list containing two elements:
    1) new string
    2) 1 if the flip was performed; 0 otherwise
    """
    new_s = s
    if new_s[i].islower():  # alpha-orbital, no swap needed;
        return new_s, 0
    else:
        j = 1
        while new_s[i + j].isupper():
            j += 1
        new_s = new_s[:i] + new_s[i + j] + new_s[i + 1:i + j] + new_s[i] + new_s[i + j + 1:]
        return new_s, 1


def place_high(s, i):
    """
    Flips two characters in the string to place the next available upper-case character at i-th position.
    E.g. place_low('aBcD',0) will return 'BacD'
    Parameters
    ----------
    s: a string containing upper and lower case characters
    i: position in the string

    Returns
    -------
    a list containing two elements:
    1) new string
    2) 1 if the flip was performed; 0 otherwise
    """
    new_s = s
    if new_s[i].isupper():  # beta-orbital, no swap needed;
        return new_s, 0
    else:
        j = 1
        while new_s[i + j].islower():
            j += 1
        new_s = new_s[:i] + new_s[i + j] + new_s[i + 1:i + j] + new_s[i] + new_s[i + j + 1:]
        return new_s, 1


def standardize_det(s):
    """
    Rearranges the orbitals in a given determinat string to the standard format,
    that is 'uLuLuL', 'uLuLuLuuu' or 'uLuLuLLLL'
    Parameters
    ----------
    s: determinant string

    Returns
    -------
    1) Determinant string in a standard format
    2) Number of pairwise flips required to obtain the standard format
    """
    new_s = s
    Nup, Ndown = 0, 0
    for i in range(len(s)):
        if s[i].isupper():
            Ndown += 1
        else:
            Nup += 1

    arr_up = 0
    arr_down = 0
    complete = False
    i = 0
    flips = 0
    while arr_down < Ndown and arr_up < Nup:
        new_s, flip = place_low(new_s, i)
        flips += flip
        arr_up += 1
        i += 1

        new_s, flip = place_high(new_s, i)
        flips += flip
        i += 1
        arr_down += 1
    return new_s, flips
