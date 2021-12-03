def delta(c1, c2):
    # Helper function used for scoring
    if c1 == c2:
        return 0
    else:
        return 1


def SP(chars, f=delta):
    total = 0
    for i in range(len(chars) - 1):
        for j in range(i + 1, len(chars)):
            total += f(chars[i], chars[j])

    return total


def SP_alignment(strings, f=delta):
    total = 0
    for character in range(len(strings[0])):
        chars = [s[character] for s in strings]
        total += SP(chars, f)
    return total

if __name__ == '__main__':
    def f(x, y):
        if x == y:
            return 2
        else:
            return -2
    print(SP(['_', 'T', '_', 'T'], f))
    print(SP_alignment(['ACG__GAGA', '_CGTTGACA', 'AC_T_GA_A', 'CCGTTCAC_'], f))
