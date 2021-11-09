from EditDistance import *
import numpy as np


def compute_indx(Sc,s2):
    indx = []
    for i in range(len(s2)):
        if i >= len(Sc) or (s2[i] == '_' and Sc[i] != '_'):
            indx.append(i)
            Sc = Sc[:i] + '_' + Sc[i:]
    indx += [len(Sc) + i for i in range(len(s2) - len(Sc))]
    return indx


def center_star(strings):
    """
    :param strings: List of strings to compute
    :return: center string, index of center string in strings array (param)
    """

    mat = np.zeros((len(strings), len(strings)))  # make k x k matrix
    for i, S in enumerate(strings):
        for j, T in enumerate(strings):
            if j <= i:  # skip the lower diagonal elements
                continue
            D = needleman(S, T)[0]
            mat[i][j] = D[len(S)][len(T)]  # get the edit distance
            mat[j][i] = mat[i][j]

    # the string that minimizes sum of distances to all other strings
    center = np.argmin(mat.sum(axis=1))
    return strings[center], center


def MSA(strings):
    Sc, center = center_star(strings)
    alignment = [s for s in strings]  # copy the input strings
    for i in range(len(strings)):
        if i == center:  # skip center string
            continue

        V, P = needleman(Sc, strings[i])
        s2, t2 = construct_alignment(P, Sc, strings[i])
        alignment[i] = t2

        add_spaces = compute_indx(Sc, s2)
        for indx in add_spaces:
            for j in range(i):
                alignment[j] = alignment[j][:indx] + '_' + alignment[j][indx:]
        Sc = s2
    alignment[center] = Sc

    return alignment


if __name__ == "__main__":
    print(center_star(['bananan', 'bann', 'aannnana', 'ananana', 'bananan']))
    msa = MSA(['bananan', 'bann', 'aannnana', 'ananana', 'bananan'])
    for s in msa:
        print(s)
