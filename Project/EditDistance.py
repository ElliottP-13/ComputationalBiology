import sys
import numpy as np
from Utils import *

# This is copied from HW2 with some tweaks to make it easier for Center_Star
# Please look at HW2 for the complete implementation


def construct_alignment(P, S, T, i=-1, j=-1, blank_s=None, blank_t=None):
    # Reconstruct the optimal alignment from the output of needleman(S,T)
    if i == -1 and j == -1:
        i, j = len(S), len(T)

    if i == 0 and j == 0:  # base case
        # print('*' * 80)  # print out the optimal alignment
        # print(out_S)
        # print(out_T)
        return blank_s, blank_t

    move = P[i][j]
    if blank_s is None or blank_t is None:
        blank_s = []
        blank_t = []

    if move == 1:  # diagonal
        i = i -1
        j = j -1
    elif move == 2:  # delete
        blank_t.append(j)
        i = i - 1
        j = j
    elif move == 3:  # insert
        blank_s.append(i)
        i = i
        j = j - 1

    blank_s, blank_t = construct_alignment(P, S, T, i, j, blank_s, blank_t)  # recur

    return blank_s, blank_t  # return the last optimal alignment


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

    for r in range(1, len(V)):  # 2nd row onwards
        for c in range(1, len(V[r])):  # 2nd col onwards

            replace = V[r-1, c-1] + f(S[r-1], T[c-1])
            delete = V[r-1, c] + f(S[r-1], '_')
            insert = V[r, c - 1] + f('_', T[c-1])

            arr = np.array([replace, delete, insert])

            V[r, c] = min(arr)
            if replace == min(arr):
                P[r][c] = 1
            if delete == min(arr):
                P[r][c] = 2
            if insert == min(arr):
                P[r][c] = 3

    return V, P


if __name__ == '__main__':
    S, T = 'ACTGGGAAA','CTGGAACA'
    V, P = needleman(S, T)
    print(f"Optimal Score: {V[-1, -1]}")
    print(V)
    blank_s, blank_t = construct_alignment(P, S, T)
    blank_s.reverse()
    blank_t.reverse()

    insert = lambda char, i, string: string[:i] + char + string[i:]
    retS = S
    retT = T
    print(type(blank_s))
    for i in blank_s:  # go backwards so don't mess up index
        retS = insert('_', i, retS)
    for i in blank_t:  # go backwards so don't mess up index
        retT = insert('_', i, retT)

    print(retS)
    print(retT)