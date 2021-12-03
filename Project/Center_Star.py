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


def center_star(strings, f=delta):
    """
    :param strings: List of strings to compute
    :return: center string, index of center string in strings array (param)
    """

    mat = np.zeros((len(strings), len(strings)))  # make k x k matrix
    for i, S in enumerate(strings):
        for j, T in enumerate(strings):
            if j <= i:  # skip the lower diagonal elements
                continue
            D = needleman(S, T, f)[0]
            mat[i][j] = D[len(S)][len(T)]  # get the edit distance
            mat[j][i] = mat[i][j]

    # the string that minimizes sum of distances to all other strings
    center = np.argmin(mat.sum(axis=1))
    return strings[center], center


def MSA(strings, f=delta):
    Sc, center = center_star(strings)
    alignment = [s for s in strings]  # copy the input strings
    for i in range(len(strings)):
        if i == center:  # skip center string
            continue

        V, P = needleman(Sc, strings[i], f)
        s2, t2 = construct_alignment(P, Sc, strings[i])
        alignment[i] = t2

        add_spaces = compute_indx(Sc, s2)
        for indx in add_spaces:
            for j in range(i):
                alignment[j] = alignment[j][:indx] + '_' + alignment[j][indx:]
        Sc = s2
    alignment[center] = Sc

    return alignment


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


def diff_generator(a, b):
    def f(c1, c2):
        # Helper function used for scoring
        if c1 == c2:
            return 0
        elif c1 == '_' or c2 == '_':
            return b
        else:
            return a
    return f


if __name__ == "__main__":
    if len(sys.argv) == 4:
        alpha = float(sys.argv[1])
        beta = float(sys.argv[2])
        file = sys.argv[3]
        seq = read_fasta(file)
        f = diff_generator(alpha, beta)
        ins = list(seq.values())

        print(f'Center (string, index): {center_star(ins, f)}')
        print("MSA:")
        msa = MSA(ins, f)
        for s in msa:
            print(s)
    else:
        print(center_star(['ACTCTCGATC', 'ACTTCGATC', 'ACTCTCTATC', 'ACTCTCTAATC']))
        print(center_star(['bananan', 'bann', 'aannnana', 'ananana', 'bananan']))
        msa = MSA(['bananan', 'bann', 'aannnana', 'ananana', 'bananan'])
        for s in msa:
            print(s)
