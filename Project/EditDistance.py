import sys
import numpy as np
from Utils import *

# This is copied from HW2 with some tweaks to make it easier for Center_Star
# Please look at HW2 for the complete implementation


def construct_alignment(P, S, T, i=-1, j=-1, out_S='', out_T=''):
    # Reconstruct the optimal alignment from the output of needleman(S,T)
    if i == -1 and j == -1:
        i, j = len(S), len(T)

    if i == 0 and j == 0:  # base case
        # print('*' * 80)  # print out the optimal alignment
        # print(out_S)
        # print(out_T)
        return out_S, out_T

    retS, retT = '', ''
    moves = P[i][j]
    for move in moves:
        if move == 0:  # is just a part of the output of P, due to its construction. So skip move
            continue
        elif move == 1:  # diagonal
            out_S = S[i - 1] + out_S
            out_T = T[j - 1] + out_T
            i = i -1
            j = j -1
        elif move == 2:  # delete
            out_S = S[i - 1] + out_S
            out_T = '_' + out_T
            i = i - 1
            j = j
        elif move == 3:  # insert
            out_S = '_' + out_S
            out_T = T[j-1] + out_T
            i = i
            j = j - 1
        retS, retT = construct_alignment(P, S, T, i, j, out_S, out_T)  # recur

    return retS, retT  # return the last optimal alignment


def needleman(S, T, f=delta):
    """
    Compute the edit distance between two strings using Needleman-Wunch Alg
    :param S: First String
    :param T: Second String
    :param f: distance function for pairwise distance. Note this implementation always minimizes distance
    :return:
    """
    V = np.zeros((len(S) + 1, len(T) + 1))  # value
    P = np.zeros((len(S) + 1, len(T) + 1))  # path

    # initialize V
    V[:, 0] = [i for i in range(len(S) + 1)]
    V[0, :] = [i for i in range(len(T) + 1)]
    P[:, 0] = 2 * np.ones(len(S) + 1)  # vertical case 2
    P[0, :] = 3 * np.ones(len(T) + 1)  # horizontal case 1
    P[0, 0] = 0  # reset 0, 0

    P = [[[x] for x in row] for row in P ]  # convert to list of lists
    for r in range(1, len(V)):  # 2nd row onwards
        for c in range(1, len(V[r])):  # 2nd col onwards

            replace = V[r-1, c-1] + f(S[r-1], T[c-1])
            delete = V[r-1, c] + f(S[r-1], '_')
            insert = V[r, c - 1] + f('_', T[c-1])

            arr = np.array([replace, delete, insert])

            V[r, c] = min(arr)
            if replace == min(arr):
                P[r][c].append(1)
            if delete == min(arr):
                P[r][c].append(2)
            if insert == min(arr):
                P[r][c].append(3)

    return V, P