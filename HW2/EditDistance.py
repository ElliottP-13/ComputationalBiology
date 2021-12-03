import sys
import numpy as np


def delta(c1, c2):
    # Helper function used for scoring
    if c1 == c2:
        return 0
    else:
        return -1


def construct_alignment(P, S, T, i=-1, j=-1, out_S='', out_T=''):
    # Reconstruct the optimal alignment from the output of needleman(S,T)
    if i == -1 and j == -1:
        i, j = len(S), len(T)

    if i == 0 and j == 0:  # base case
        print('*' * 80)  # print out the optimal alignment
        print(out_S)
        print(out_T)
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


def needleman(S, T):
    V = np.zeros((len(S) + 1, len(T) + 1))  # value
    P = np.zeros((len(S) + 1, len(T) + 1))  # path

    # initialize V
    V[:, 0] = [-i for i in range(len(S) + 1)]
    V[0, :] = [-i for i in range(len(T) + 1)]
    P[:, 0] = 2 * np.ones(len(S) + 1)  # vertical case 2
    P[0, :] = 3 * np.ones(len(T) + 1)  # horizontal case 1
    P[0, 0] = 0  # reset 0, 0

    P = [[[x] for x in row] for row in P ]  # convert to list of lists
    for r in range(1, len(V)):  # 2nd row onwards
        for c in range(1, len(V[r])):  # 2nd col onwards

            replace = V[r-1, c-1] + delta(S[r-1], T[c-1])
            delete = V[r-1, c] + delta(S[r-1], '_')
            insert = V[r, c - 1] + delta('_', T[c-1])

            arr = np.array([replace, delete, insert])

            V[r, c] = max(arr)
            if replace == max(arr):
                P[r][c].append(1)
            if delete == max(arr):
                P[r][c].append(2)
            if insert == max(arr):
                P[r][c].append(3)

    return V, P


def read_fasta(filepath):
    """
    Returns a dict key=sequence_name, value=sequence from the fasta file
    :param filepath:
    :return:
    """
    file = open(filepath)
    print(f"Opening: {filepath}")
    started = False
    header = ''
    text = ''
    sequences = {}
    for line in file:  # reads in fasta file. Allows for sequence to continue on multiple lines
        if line.startswith('>'):
            if started:
                sequences[header] = text.replace('\n', '')
                text = ''
            started = True
            header = line.replace('>', '').strip()
        else:
            text += line.strip()
    sequences[header] = text.replace('\n', '')
    return sequences


if __name__ == '__main__':
    if len(sys.argv) == 2:
        # run fasta file
        # expect arguments z_algorithm.py <path to file>
        filepath = sys.argv[1]
        seq = read_fasta(filepath)
        vals = list(seq.values())
        keys = list(seq.keys())
        print(f"Computing alignment between S={keys[0]} and T={keys[1]}")
        V, P = needleman(vals[0], vals[1])
        print(f"Optimal Score: {V[-1,-1]}")
        print(V)
        s, t = construct_alignment(P, vals[0], vals[1])
    else:
        V, P = needleman('ACTCTCGATC', 'ACTTCGATC')
        print(f"Optimal Score: {V[-1, -1]}")
        print(V)
        s, t = construct_alignment(P, 'ACTCTCGATC', 'ACTTCGATC')
        pass
