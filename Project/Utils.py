def delta(c1, c2):
    # Helper function used for scoring
    if c1 == c2:
        return 0
    else:
        return 1


def SP(chars):
    total = 0
    for i in range(len(chars) - 1):
        for j in range(i + 1, len(chars)):
            total += delta(chars[i], chars[j])

    return total


def SP_alignment(strings):
    total = 0
    for character in range(len(strings)):
        chars = [s[character] for s in strings]
        total += SP(chars)
    return total
