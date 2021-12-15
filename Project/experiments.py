import Center_Star as cs
import ClustalW as cw
import ExactDP as ed
from EditDistance import *
from Utils import *
import random
import time


def msa_to_cell(msa):
    o = '\"'
    o += ', '.join(msa)
    return o + '\"'


def run_experiment(strings, file, run_exact, func, threshold=900):
    if run_exact:
        st = time.time()
        exact_score, (V, P) = ed.generalized_needleman(strings, f=func, compute_path=True)
        exact_time = time.time() - st

        exact_msa = ed.construct_alignment(P, strings)
        file.write(f'{exact_time}, {exact_score}, {msa_to_cell(exact_msa)}, ')
        print(f'{exact_time}, {exact_score}, {msa_to_cell(exact_msa)}, ')

        run_exact = exact_time < threshold
    else:
        print(f'Skipping exact')
        file.write(' , , , ')

    center_star = cs.CenterStar(func)

    st = time.time()
    cs_msa = center_star.MSA(strings)
    cs_time = time.time() - st

    cs_score = SP_alignment(cs_msa, func)

    file.write(f'{cs_time}, {cs_score}, {msa_to_cell(cs_msa)}, ')
    print(f'{cs_time}, {cs_score}, {msa_to_cell(cs_msa)}, ')

    st = time.time()
    cw_msa = cw.clustalW(strings, draw=False, f=func)
    cw_time = time.time() - st

    cw_score = SP_alignment(cw_msa, func)

    file.write(f'{cw_time}, {cw_score}, {msa_to_cell(cw_msa)}\n')
    print(f'{cw_time}, {cw_score}, {msa_to_cell(cw_msa)}\n')
    file.flush()

    return run_exact


def run_rand_strings_experiments(out_f = './results/num_strings.csv', string_len = 50, seed=0, func=delta):
    alphabet = ['A', 'C', 'T', 'G']
    random.seed(seed)

    # Once a trial exceeds this time, we will no longer run the exact algorithm
    # Note, we do not terminate the trial, but let it finish, so the last trial will exceed this time

    file = open(out_f, 'a')
    file.write('k, exact time, exact score, exact align, star time, star score, star align, '
               'clustal time, clustal score, clustal align\n')

    run_exact = True

    for k in range(3, 50):
        print(f'K = {k}')
        strings = [random.choices(alphabet, k=string_len) for _ in range(k)]  # generate k random words

        file.write(f'{k}, ')

        run_exact = run_experiment(strings, file, run_exact, func)




def run_center_strings_experiments(out_f = './results/num_strings.csv', string_len = 50, seed=0, func=delta):
    alphabet = ['A', 'C', 'T', 'G']
    random.seed(seed)

    run_exact = True
    # Once a trial exceeds this time, we will no longer run the exact algorithm
    # Note, we do not terminate the trial, but let it finish, so the last trial will exceed this time

    file = open(out_f, 'a')
    file.write('k, exact time, exact score, exact align, star time, star score, star align, '
               'clustal time, clustal score, clustal align\n')

    for k in range(3, 50):
        print(f'K = {k}')
        base = random.choices(alphabet, k=string_len)  # generate base string
        strings = [base]
        for i in range(1, k):
            s = base
            for mut in range(random.randint(5, string_len // 3)):  # random number of mutations
                j = random.randrange(0, string_len)
                s = s[:j] + random.choice(alphabet) + s[j+1:]
                strings.append(s)

        file.write(f'{k}, ')

        run_exact = run_experiment(strings, file, run_exact, func)


if __name__ == "__main__":
    run_rand_strings_experiments(out_f='./results/random_strings_len50.csv', string_len=50)
    run_center_strings_experiments(out_f='./results/center_strings_len50.csv', string_len=50)
    pass
