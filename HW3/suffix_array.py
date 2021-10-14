
def suffix_array(S):
    suffix = [(i, S[i:]) for i in range(len(S))]  # generate list of tuples
    suffix.sort(key=lambda tup: tup[1])  # sort by suffix's
    return [fix[0] for fix in suffix]

if __name__ == '__main__':
    print(suffix_array("banana$"))
    pass