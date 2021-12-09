import collections
import networkx as nx
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


def make_delta_fn(S, T, f=delta):
    """
    Makes PSP Function
    PSP(i,j) = sum_{x,y} g_x(i) g_y(j) delta(x, y)
    where g_c(i) indicates frequency of char c in column i
    :param S:
    :param T:
    :param f:
    :return:
    """
    # get list of dicts. list[i] = dict{char: frequency in column}
    s_char = [collections.Counter([st[i] for st in S]) for i in range(len(S[0]))]
    s_char.append({'_': 1})
    t_char = [collections.Counter([tt[i] for tt in T]) for i in range(len(T[0]))]
    t_char.append({'_': 1})

    def g(i, j):
        total = 0
        for x in s_char[i].keys():
            for y in t_char[j].keys():
                total += s_char[i][x] * t_char[j][y] * f(x, y)
        return total

    return g


def merge_on_tree(strings, T, node_str=None):
    """
    Merge the set of strings using the guide tree
    :param strings: set of strings for MSA
    :param T: guide tree
    :param node_str: recursive dict mapping aligned sequence at each node
    :return: MSA alignment
    """
    if node_str is None:
        node_str = {}
        for i in range(len(strings)):
            node_str[i + 1] = [strings[i]]  # stores in mapping

    leaves = [x for x in T.nodes() if T.degree(x) == 1]  # get all the leaf nodes
    one_parent = {}  # temp holding list
    parents = []
    for child in leaves:
        edge = list(T.edges(child))[0]
        if edge[1] not in one_parent:
            one_parent[edge[1]] = child
        else:  # add it to the list if there are 2 children with this parent
            c2 = one_parent.pop(edge[1])  # get parent and 2nd child
            parents.append((edge[1], child, c2))  # (parent, child 1 , child 2)

    if len(parents) == 0:  # recursive case
        return node_str[0]

    for p, c1, c2 in parents:  # for ones we can merge
        A = node_str[c1]  # set of strings S
        B = node_str[c2]  # set of strings T

        f = make_delta_fn(A, B)  # make PSP function
        V, P = needleman(A[0], B[0], f, True)  # assume everything in A, B is same length respectively

        blank_s, blank_t = construct_alignment(P)
        blank_s.reverse()
        blank_t.reverse()  # go backwards so don't mess up index

        insert = lambda char, i, string: string[:i] + char + string[i:]
        retS = A.copy()
        retT = B.copy()

        # add the blanks
        for k in range(len(retS)):
            for i in blank_s:
                retS[k] = insert('_', i, retS[k])
        for k in range(len(retT)):
            for i in blank_t:
                retT[k] = insert('_', i, retT[k])

        node_str[p] = retS + retT

        # remove edges
        T.remove_edge(c1, p)
        T.remove_edge(c2, p)

        # remove nodes
        T.remove_node(c1)
        T.remove_node(c2)

    return merge_on_tree(strings, T, node_str)


def clustalW(strings, draw=False):
    distance_matrix = compute_pairwise_distances(strings)
    T = neighbor_join(distance_matrix)

    if draw:
        pos = nx.drawing.planar_layout(T)
        nx.draw(T, pos=pos, with_labels=True)
        plt.show()

    MSA = merge_on_tree(strings, T)
    return MSA


if __name__ == '__main__':
    msa = clustalW(['ACTCTCGATC', 'ACTTCGATC', 'ACTCTCTATC', 'ACTCTCTAATC'])
    for s in msa:
        print(s)
    print(f' ANSWER: {SP_alignment(msa)}')
    print(f" ONLINE ANSER: {SP_alignment(['ACTCTCTATC_', 'ACTCTCTAATC', 'ACTCTCGATC_', 'ACT_TCGATC_'])}")
    """
    'ACTCTCTATC_', 'ACTCTCTAATC', 'ACTCTCGATC_', 'ACT_TCGATC_'
    """

