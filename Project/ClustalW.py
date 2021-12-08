import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from EditDistance import *


def compute_pairwise_distances(strings, f=delta):
    """
    Compute all pairwise distances using needleman wunch alg
    :param strings: Strings to compare
    :param f: score function
    :return: distance matrix
    """
    mat = np.zeros((len(strings), len(strings)))
    for i in range(len(strings)):
        for j in range(i + 1, len(strings)):
            mat[i][j] = needleman(strings[i], strings[j], f)[0][-1][-1]
            mat[j][i] = mat[i][j]

    return mat


def neighbor_join(D, T=None, m=None):
    """
    Implement the neighbor join algorithm
    :param D: distance matrix
    :return: Tree of merged results
    """
    if T is None:
        T = nx.star_graph(len(D))  # makes 1 leaf for every string. Node 0 is the center
        m = [i+1 for i in range(len(D))]  # map index in D to node in T

    if len(D) == 2:
        return T

    # Based on Wikipedia algorithm description
    # Step 1 make Q matrix
    Q = np.zeros(D.shape)
    for i in range(len(Q)):
        for j in range(len(Q[i])):
            s1 = np.sum(D[i][:])
            s2 = np.sum(D[j][:])
            Q[i][j] = (len(D) - 2) * D[i][j] - s1 - s2
    print(Q)
    np.fill_diagonal(Q, np.infty)
    # Step 2 find pair with smallest value to join into u
    idx = np.argmin(Q)  # index into flat array
    i, j = np.unravel_index(idx, Q.shape)  # convert into i, j pair
    mi, mj = m[i], m[j]  # map onto the nodes

    # Update T
    # nx.draw(T, with_labels=True)
    # plt.show()
    T.remove_edge(mi, 0)

    # nx.draw(T, with_labels=True)
    # plt.show()
    T.remove_edge(mj, 0)

    mu = T.number_of_nodes()

    # compute edge weight i, u
    s1 = np.sum(D[i][:])
    s2 = np.sum(D[j][:])
    w = (1/2) * (len(D) - 2) * D[i][j] + (1 / (2 * (len(D)-2)))*(s1 - s2)
    T.add_node(mu)
    T.add_edge(mi, mu, weight=w)
    T.add_edge(mj, mu, weight=(D[i][j] - w))
    T.add_edge(mu, 0)

    u = len(Q)  # last index
    m.append(mu)  # store the mapping

    # Step 3 compute distance to new node u
    Dnew = np.column_stack((D, np.zeros(len(D))))  # add blank col
    Dnew = np.vstack((Dnew, np.zeros(len(Dnew[0]))))  # add blank row

    # fill in new spots
    for k in range(len(D)):
        d = 1/2 * (D[i][k] + D[j][k] - D[i][j])
        Dnew[u][k] = d
        Dnew[k][u] = d

    # Delete rows and cols i, j
    Dnew = np.delete(np.delete(Dnew, min(i,j), 1), min(i,j), 0)  # remove one
    Dnew = np.delete(np.delete(Dnew, max(i,j)-1, 1), max(i,j)-1, 0)  # remove other

    # shift the mappings
    for a in range(min(i, j) + 1, len(m)):  # shift left by 1
        m[a - 1] = m[a]
    for a in range(max(i, j), len(m)):  # shift left by 1
        m[a - 1] = m[a]

    return neighbor_join(Dnew, T, m[:-2])


def clustalW(strings):
    distance_matrix = compute_pairwise_distances(strings)
    T = neighbor_join(distance_matrix)
    pos = nx.drawing.planar_layout(T)
    nx.draw(T, pos=pos, with_labels=True)
    plt.show()


if __name__ == '__main__':
    clustalW(['actct','agcat','agct','acttg','ctct'])
