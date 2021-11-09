from EditDistance import *
import numpy as np


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


if __name__ == "__main__":
    print(center_star(['bananan', 'bann', 'aannnana', 'ananana', 'bananan']))
    pass
