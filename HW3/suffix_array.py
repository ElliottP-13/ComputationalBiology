
def suffix_array(S):
    suffix = [(i, S[i:]) for i in range(len(S))]  # generate list of tuples
    suffix.sort(key=lambda tup: tup[1])  # sort by suffix's
    return [fix[0] for fix in suffix]


def pattern_search(pattern, S):
    suffix = suffix_array(S)
    l, r = 0, len(S)-1
    found = False
    while l <= r:
        m = int((l + r) / 2)
        st = S[suffix[m]:]  # string version of the suffix
        if pattern == st[:len(pattern)]:
            print(f'Found match at {suffix[m]}')
            found = True
            break
        elif pattern > st:
            l = m + 1
        else:
            r = m - 1

    if not found:
        print("Pattern not found in string")

if __name__ == '__main__':
    print(suffix_array("banana$"))
    pattern_search('ba', 'banana$')