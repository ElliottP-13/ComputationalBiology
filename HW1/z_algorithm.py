"""
Author - Elliott Pryor
Created on - 2 September 2021
"""
import sys


def print_output(z):
    if len(z) > 0:
        z = [str(i) for i in z]
        output = ', '.join(z)
        print(f"Found matches {output}")
    else:
        print("No match found")

def find_next_match(s1, s2, S):
    """
    Finds the next match in S starting from s1 and s2
    Compares S[s1, n] to S[s2, n]
    We don't pass in slices because this is faster
    :param s1: Starting of first string (base)
    :param s2: Start of second string (prefix to match)
    :param S: String to search
    :return: q=s2 index of failed match in S, zk length of match
    """
    n = len(S)
    zk = 0
    while (s1 < n and s2 < n) and (S[s1] == S[s2]):  # while they match and in bounds
        s1 += 1
        s2 += 1
        zk += 1
    return s2, zk


def z_alg(S):
    """
    Run the z prefix matching algorith
    :param S: String input
    :return: Z score for each index (note first one is always 0)
    """
    l, r = 0, 0  # initialize
    z = [0 for _ in S]  # init z scores
    for k in range(1, len(S)):
        if k > r:  # match them
            _, zk = find_next_match(0, k, S)
            if zk > 0:
                r = k + zk - 1
                l = k
                z[k] = zk
        else:
            kp = k - l
            b = r - k + 1  # ||B|| from alg descrip
            if z[kp] < b:
                z[k] = z[kp]
            else:
                q, zk = find_next_match(b, r + 1, S)
                z[k] = q - k
                l = k
                r = q - 1
    return z


def z_pattern_match(pattern, text, special_char=None):
    """
    Finds all occurrences of the pattern in the text
    :param pattern: pattern to search for
    :param text: text to search in
    :param special_char: special division characters. Can specify a char, or list of chars.
                         Defaults are $, %, ^, ~, `, #
    :return: Index in text of start of pattern match
    """
    if special_char is None:
        special_char = ['$', '%', '^', '@', '~', '`', '#']  # our default options for special chars

    sep_char = ''
    if not isinstance(special_char, str):  # if it isn't just a string
        for c in special_char:
            if c not in text and c not in pattern:  # 2 * O(n) each time
                sep_char = c
                break
    S = pattern + sep_char + text  # put special char in the middle
    z = z_alg(S)
    m = len(pattern)
    matches = []
    for j in range(1, len(z)):
        if z[j] == m:
            matches.append(j - (m + 1))
    return matches


def analyze_fasta(pattern, filepath):
    file = open(filepath)
    print(f"Opening: {filepath}")
    print(f"Matching pattern: {pattern}")
    started = False
    text = ''
    for line in file:
        if line.startswith('>'):
            if started:
                z = z_pattern_match(pattern, text.replace('\n', ''))
                print_output(z)
            started = True
            print(f"Analyzing Sequence: {line.replace('>', '').strip()}")
        else:
            text += line.strip()
    z = z_pattern_match(pattern, text.replace('\n', ''))
    print_output(z)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        # expect arguments z_algorithm.py <pattern> <path to file>
        pattern = sys.argv[1]
        filepath = sys.argv[2]
        analyze_fasta(pattern, filepath)
    else:
        print(z_pattern_match('aa', 'aaaaaaaaaaaaa'))

