import numpy as np
from EditDistance import *


def SP(chars):
    total = 0
    for i in range(len(chars) - 1):
        for j in range(i + 1, len(chars)):
            total += delta(chars[i], chars[j])

    assert sum([delta(chars[i], chars[j]) for i in range(len(chars) - 1) for j in range(i+1, len(chars))]) == total
    return total


def generalized_needleman(strings, f, compute_path=True):

    strings = ['_' + s for s in strings]  # prepend blank so we can use 1 indexing

    V = np.zeros([len(s) for s in strings])  # value
    if compute_path:
        P = np.zeros([len(s) for s in strings])  # value

    k = len(strings)

    # for each possible b combination

    for idx in np.ndindex(V.shape):
        if sum(idx) == 0:  # skip V[0]
            continue

        print(idx)
        vopt = float('inf')
        iopt = None
        for i in range(1, 2 ** k):
            b = format(i, 'b').zfill(k)
            b = np.array([int(c) for c in b])
            jdx = idx - b

            if np.any(jdx < 0):  # recurrence relation i_j = 0 --> b_j = 0
                continue

            jdx = tuple(jdx)  # to be able to index into V or P

            v = V[jdx] + SP([strings[m][jdx[m]] for m in range(len(strings))])
            if v < vopt:
                vopt = v
                iopt = i

        V[idx] = vopt
        if compute_path:
            P[idx] = iopt
    return V[tuple(x - 1 for x in V.shape)]


if __name__ == "__main__":
    s = ['ACTCTCGATC', 'ACTTCGATC', 'ACTCTCTATC', 'ACTCTCTAATC']
    q = np.array(s)

    print(f' ANSWER: {generalized_needleman(s, delta)}')
    pass
