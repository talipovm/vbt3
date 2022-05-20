import numpy
import scipy.linalg
from vbt3 import FixedPsi

# manual control of the parallelization feature
PARALLEL = False

RAY_IMPORTED = False
if PARALLEL:
    # attempt to import ray
    try:
        import ray
        RAY_IMPORTED = True
        my_decor = ray.remote
    except ImportError as e:
        pass

if not (PARALLEL and RAY_IMPORTED):
    # No parallelization; use a dummy decorator
    def my_decor(func):
        return func


def repair_connections(coupled, z):
    # Repair the coupled hash
    # Problem: some elements of z are 0 because of the numerical fluctuations
    # Solution: recover these elements using identities:
    # For a < b < c:
    # z[a,b] = v[b] / v[a]
    # z[a,c] = v[c] / v[a]
    # z[b,c] = v[c] / v[b]
    # Thus,
    # z[a,b] = z[a,c] / z[b,c]
    #
    # E.g., z[0,20] = 0 but z[0,19]=1 and z[19,20]=3.49
    # z[0,20] can be recovered as = z[0,19] / z[20,19]
    repaired = {}
    for k in coupled.keys():
        repaired[k] = {}
        for m in coupled[k].keys():
            repaired[k][m] = coupled[k][m]
            if coupled[k][m] == 0.:
                # find a non-zero valued index to use for repair, z[k][m] = z[k][nz] / z[m][nz] if m < nz
                # if m > nz, z[k][m] = z[k][nz] * z[nz][m]
                fixed = False
                for nz in coupled[k].keys():
                    if z[k, nz] != 0.0:
                        if m < nz and z[m,nz] != 0.:
                            repaired[k][m] = z[k, nz] / z[m, nz]
                            fixed = True
                        elif m > nz and z[nz, m] != 0.:
                            repaired[k][m] = z[k, nz] * z[nz, m]
                            fixed = True
                        if fixed:
                            break
    return repaired


@my_decor
def single_trial(mH, mS, h, s, precision=12):
    N = mH.shape[0]
    # get the numeric eigenvectors
    H = numpy.array(mH.subs({'h': h, 's': s})).astype(numpy.float64)
    S = numpy.array(mS.subs({'s': s})).astype(numpy.float64)
    vals, vecs = scipy.linalg.eig(H, S)
    v = vecs[:, vals.argmin(0)]

    # compute the coefficient ratios
    fixed_coefs = numpy.zeros((N, N))
    for i in range(N):
        for j in range(i + 1, N):
            if v[i] == 0.0:
                fixed_coefs[i, j] = numpy.Inf
            else:
                fixed_coefs[i, j] = numpy.around(v[j] / v[i], precision)
    return fixed_coefs


def get_coupled(mS, mH, N_tries=10, precision=12, ranges={'h': (-1.0, 0.0), 's': (0.0, 1.0)}):
    """
    Get the dictionary showing the components of the lowest energy wave vector
    that have constant ratio independent of the numerical values h, s
    """
    N = mH.shape[0]

    hv = numpy.random.uniform(low=ranges['h'][0], high=ranges['h'][1], size=N_tries)
    sv = numpy.random.uniform(low=ranges['s'][0], high=ranges['s'][1], size=N_tries)

    fcs = [None, ] * N_tries

    if RAY_IMPORTED:
        # Parallel
        ray.init()
        ids = [None, ] * N_tries
        for trial in range(N_tries):
            ids[trial] = single_trial.remote(mH, mS, hv[trial], sv[trial], precision=precision)
        # Block until the tasks are done and get the results.
        fcs = ray.get(ids)
        ray.shutdown()
    else:
        # Single core
        for trial in range(N_tries):
            fcs[trial] = single_trial(mH, mS, hv[trial], sv[trial], precision=precision)

    # check if the coefficient ratios differ from the first value
    b = numpy.ones((N, N))
    for trial in range(1, N_tries):
        b *= (fcs[0] == fcs[trial])
    z = numpy.triu(b * fcs[0], k=1)

    # Convert the matrix in the dictionary form
    coupled = {}
    for i in range(z.shape[0]):
        # i shows the row
        for j in range(i + 1, z.shape[0]):
            if z[i, j] == 0:
                continue
            # j is the index for a nonzero coupled i,j
            found = False
            for k, v in coupled.items():
                # k is the reference FixedPsi, it has weight 1.0 by definition
                # v is a dict of the combined FixedPsi index/coef pairs
                if i in v:
                    found = True
                    if j not in v:
                        coupled[k][j] = z[k, j]
                    break
            if not found:
                coupled[i] = {i: 1.0, j: z[i, j]}

    import pickle
    # file = open('z.pickle', 'wb')
    # pickle.dump(z, file)
    # file.close()
    # file = open('coupled.pickle', 'wb')
    # pickle.dump(coupled, file)
    # file.close()

    return repair_connections(coupled, z)



def get_combined(P, indices, coefs=None):
    """
    Contract an array of FixedPsi objects
    :param indices: indices of the FixedPsi objects to be combined together
    :param coefs: coefficients of the FixedPsi objectsA contraction scheme. Assumed to be 1.0 if omitted
    :return: A new array of FixedPsi objects with combined determinants
    """
    if coefs is None:
        coefs = numpy.ones_like(indices)

    P_new = []
    P_contracted = FixedPsi()

    # create a combined FixedPsi
    for i in range(len(indices)):
        P_contracted.add_FixedPsi(P[indices[i]], coefs[i])
    P_new.append(P_contracted)

    # copy the remaining FixedPsi
    for i in range(len(P)):
        if i not in indices:
            P_new.append(P[i])

    return P_new


def get_combined_from_dict(P, d):
    """
    Contract a list of FixedPsi objects
    :param P: original list of FixedPsi
    :param d: dictionary of dictionaries; each sub-dict has FixedPsi indices/coefs as key/value pairs
    :return: A new list of FixedPsi objects with combined determinants
    """

    P_new = []
    c_indices = []

    # create the combined FixedPsi
    for c in d.values():
        P_contracted = FixedPsi()
        for index, coef in c.items():
            P_contracted.add_FixedPsi(P[index], coef)
            c_indices.append(index)
        P_new.append(P_contracted)

    # copy the remaining FixedPsi
    for i in range(len(P)):
        if i not in c_indices:
            P_new.append(P[i])

    return P_new


def validate_solution(expression, mH, mS, N_tries=10, precision=12):
    """
    Validate that formula for the lowest energy matches the numerical values
    :param expression: Sympy expression for the formula
    :param mH: symbolic Hamiltonian matrix
    :param mS: symbolic overlap matrix
    :param N_tries: number of trials
    :param precision: 10^-precision is the matching threshold
    :return: Boolean value. True if the numerical and analytical solutions produce the same result
    """
    hv = numpy.random.uniform(low=-1.0, high=0.0, size=N_tries)
    sv = numpy.random.uniform(low=0.0, high=1.0, size=N_tries)

    for trial in range(N_tries):
        # Numerical result
        H = numpy.array(mH.subs({'h': hv[trial], 's': sv[trial]})).astype(numpy.float64)
        S = numpy.array(mS.subs({'s': sv[trial]})).astype(numpy.float64)
        vals, vecs = scipy.linalg.eig(H, S)
        num_solution = min(vals)

        # Result from the formula
        formula_solution = expression.evalf(subs={'h': hv[trial], 's': sv[trial]})

        if abs(num_solution - formula_solution) >= precision:
            return False
    return True
