import numpy as np
from EditDistance import *
from Utils import *


def generalized_needleman(strings, f=delta, compute_path=True):

    strings = ['_' + s for s in strings]  # prepend blank so we can use 1 indexing

    V = np.zeros([len(s) for s in strings])  # value
    if compute_path:
        P = np.zeros([len(s) for s in strings], dtype=np.int32)  # value

    k = len(strings)

    for idx in np.ndindex(V.shape):
        if sum(idx) == 0:  # skip V[0]
            continue

        # print(idx)
        vopt = float('inf')
        iopt = None
        for i in range(1, 2 ** k):
            b = format(i, 'b').zfill(k)
            b = np.array([int(c) for c in b])
            jdx = idx - b

            if np.any(jdx < 0):  # recurrence relation i_j = 0 --> b_j = 0
                continue  # we can skip if we hit this case (will be covered by other b)

            jdx = tuple(jdx)  # to be able to index into V or P
            pdx = np.multiply(idx, b)  # element-wise multiplication

            v = V[jdx] + SP([strings[m][pdx[m]] for m in range(len(strings))], f=f)
            if v < vopt:
                vopt = v
                iopt = i

        V[idx] = vopt
        if compute_path:
            P[idx] = iopt  # less memory to store i than b

    if not compute_path:
        return V[tuple(x - 1 for x in V.shape)]
    else:
        return V[tuple(x - 1 for x in V.shape)], (V, P)


def construct_alignment(P, strings):
    k = len(strings)
    ret_strings = ['' for _ in strings]  # copy the list of strings

    index = np.prod([s for s in P.shape]) - 1  # get last index in flat array
    index = np.unravel_index(index, P.shape)

    while sum(index) > 0:
        p = P[index]
        bvec = format(p, 'b').zfill(k)
        bvec = np.array([int(c) for c in bvec])  # get the b vector that we used

        for i, b in enumerate(bvec):
            if b == 0:  # add blank
                ret_strings[i] = '_' + ret_strings[i]
            else:
                ret_strings[i] = strings[i][index[i] - 1] + ret_strings[i]

        index = tuple(index - bvec)
    return ret_strings


if __name__ == "__main__":
    s = ['ACTCTCGATC', 'ACTTCGATC', 'ACTCTCTATC', 'ACTCTCTAATC']
    q = np.array(s)

    score, (V, P) = generalized_needleman(s, compute_path=True)

    print(f' ANSWER: {score}')

    msa = construct_alignment(P, s)
    for s in msa:
        print(s)
    print(SP_alignment(msa))

    pass
