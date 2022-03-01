from prime import prime_range, isprime
from utils import smooth_factor
from factor import factor

from math import sqrt, log, exp, gcd
from random import Random
from linalg import matrix_copy
from dlp import pollard_rho_dlp


required_relations_ratio = 2

rand = Random()


def row_reduce_mod(a, row, col, m):
    h = len(a)
    w = len(a[row])
    for i in range(h):
        if i != row:
            for j in range(w):
                a[i][j] = (a[i][j] - a[row][j] * a[i][col]) % m


def make_pivot_mod(a, index, m):
    inv = pow(a[index], -1, m)
    return list(map(lambda n: (n * inv) % m, a))


def solve_matrix(matrix, n):
    print(f"isprime {n}: {isprime(n)}")

    pivot_row = 0  # first pivot belongs in first row

    w = len(matrix[0]) - 1
    h = len(matrix)

    array = matrix_copy(matrix)

    for j in range(w):

        # start at looking for pivot after previous pivot row
        for i in range(pivot_row, h):

            # if non-zero element, this row can become pivot row
            if array[i][j]:

                if array[i][j] != 1:
                    # Make j'th element the pivot, reducing rest of row as well
                    array[i] = make_pivot_mod(array[i], j, n)
                    assert array[i][j] == 1

                if i > pivot_row:  # if pivot row not already in correct position, swap
                    array[i], array[pivot_row] = array[pivot_row][:], array[i][:]

                row_reduce_mod(array, pivot_row, j, n)      # row reduce everything else
                pivot_row += 1
                break

    print(f"Array: \n{array}")

    if any(array[0][j] for j in range(1, w)):
        print(f"System was not fully solved")
        exit(0)
    else:
        solutions = []
        for i in range(h):
            if not any(array[i][j] for j in range(w)):
                return solutions
            solutions.append(array[i][w])

        return solutions


def precomp(g, n):
    global required_relations_ratio

    B = int(exp(sqrt(log(n) * log(log(n)))))

    factor_base = prime_range(B)

    required_relations = int(len(factor_base) * required_relations_ratio)

    smooth_matrix = []
    relations_found = 0

    order = n - 1
    while relations_found < required_relations:
        k = rand.randrange(n - 1)
        if (powers := smooth_factor(pow(g, k, n), factor_base)) is not None:
            smooth_matrix.append([e % order for e in powers] + [k])
            relations_found += 1

    print(f"Matrix:\n{smooth_matrix}")
    if isprime(order):
        logs = solve_matrix(smooth_matrix, order)
        for p, l in zip(factor_base, logs):
            print(f"log_g({p}) = {l % order}/{-l % order}; rho: {pollard_rho_dlp(g, p, n) % order}")
            # print(f"{p} = {pow(g, l % order, n)} or {pow(g, -l % order, n)}")
    else:
        factors = factor(order)


if __name__ == "__main__":
    precomp(4, 74)
