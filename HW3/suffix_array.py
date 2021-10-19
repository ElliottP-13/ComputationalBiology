import sys


def suffix_array(S):
    suffix = [(i, S[i:]) for i in range(len(S))]  # generate list of tuples
    suffix.sort(key=lambda tup: tup[1])  # sort by suffixes
    return [fix[0] for fix in suffix]


def pattern_search(pattern, S):
    suffix = suffix_array(S)
    l, r = 0, len(S)-1
    found = False
    while l <= r:
        m = int((l + r) / 2)
        st = S[suffix[m]:]  # string version of the suffix
        if pattern == st[:len(pattern)]:
            print(f'Found a match at {suffix[m]}')
            found = True
            break
        elif pattern > st:
            l = m + 1
        else:
            r = m - 1

    if not found:
        print("Pattern not found in string")


if __name__ == '__main__':
    pattern = None

    if len(sys.argv) == 2:
        input_str = sys.argv[1]
    elif len(sys.argv) == 3:
        input_str = sys.argv[1]
        pattern = sys.argv[2]
    else:
        input_str = input("Input the string to compute a suffix array for: ")

    if pattern is None:
        pattern = input("Input the pattern to search for: ")

    print(f"Computing suffix array for: {input_str}")
    print(suffix_array(input_str))
    print(f"Finding the pattern: {pattern}")
    pattern_search(pattern, input_str)