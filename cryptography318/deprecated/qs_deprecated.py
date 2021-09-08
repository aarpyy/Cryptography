
@deprecated
def __factor_perfect_square(n, B=7):
    """Attempts a similar attack to quadratic sieve, finding B-smooth perfect squares in an attempt
    to find a multiple of n. Function written for learning purposes, if trying to factor integer, using
    FactorInt."""

    from itertools import combinations_with_replacement as _all

    m = prime_pi(B)
    b_smooth_nums = []
    squared_nums = {}
    a = isqrt(n) - 1
    while len(b_smooth_nums) < m:
        ci = pow(a, 2, n)
        if b_smooth(ci, B):
            b_smooth_nums.append(ci)
            squared_nums[ci] = a
        a += 1

    ci_factors = {}
    factor_base = primes_lt(B)
    mat = []
    for num in b_smooth_nums:
        ci = num
        exp = []
        for p in factor_base:
            count = 0
            while num % p == 0:
                num //= p
                count += 1
            exp.append(count)
        ci_factors[ci] = exp
        mat.append(exp)

    for c in range(2, len(mat)):
        possible = list(map(list, list(_all(mat, c))))
        for comb in possible:
            acc = comb[0][:]
            perfect_2 = True

            for i in range(1, len(comb)):
                for j in range(len(acc)):
                    acc[j] += comb[i][j]
            for j in acc:
                if j % 2 != 0:
                    perfect_2 = False

            if not perfect_2:
                continue
            a, b = 1, 1
            for num in b_smooth_nums:
                choice = ci_factors[num]
                if choice in comb:
                    a *= squared_nums[num]
            for k in range(len(acc)):
                b *= pow(factor_base[k], acc[k] // 2)

            p, q = gcd(a - b, n), gcd(a + b, n)
            if 1 < p < n or 1 < q < n:
                return _factor_with_known(p, q, n)

    return False


@deprecated
def gaussian_elimination_mod(matrix):
    """Performs Gaussian elimination mod 2 over a matrix.

    Credit to https://github.com/mikolajsawicki/quadratic-sieve/tree/main/quadratic_sieve for technique of
    performing gaussian elimination on the transpose of the binary matrix."""

    m = matrix.copy().transpose()

    pivot = 0
    for j in range(len(m[0])):
        for i in range(pivot, r := len(m)):
            if pivot == r:
                break

            if m[i, j] == 1:

                # if row with pivot not in expected pivot row, swap them
                if i > pivot:
                    temp = m[i].copy()
                    m[i] = m[pivot].copy()
                    m[pivot] = temp

                # row reduce down entire matrix
                for k in range(len(m)):
                    if k == pivot:
                        continue
                    if m[k, j] == 1:
                        m[k] = (m[k] + m[pivot]) % 2
                pivot += 1

    m.remove_null_row()

    return m


@deprecated
def kernel(matrix):
    """Finds the basis of the kernel of a binary matrix. Intended as helper function for quadratic sieve algorithm.
    Assumes input matrix is in rref, has no null-rows, and is over the field Z-2 (field of integers mod 2).

    :param matrix: list[list[int]], usually Matrix or numpy.array object
    :return: dict w/
    keys: columns of free variables
    values: columns of dependent variables"""

    # finds pivots, stores in dictionary key: row index, value: col index
    pivots = {}
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if matrix[i][j] == 1:
                pivots[i] = j
                break

    # gets list of pivot columns into a set
    columns_w_pivots = set(pivots.values())
    kernel_basis = {}
    null_cols = set()
    for j in range(width := len(matrix[0])):
        if j in columns_w_pivots:
            continue

        non_zero = set()

        # iterate down matrix, if find non-zero, then this column will have a basis vector in kernel
        for i in range(len(matrix)):
            if matrix[i][j] == 1:
                non_zero.add(i)

        if len(non_zero) >= 1:
            kernel_basis[j] = [[0] * width]
            for c in non_zero:
                kernel_basis[j][0][pivots[c]] = 1

            # add current column index to basis vector as well
            kernel_basis[j][0][j] = 1

        # if no non-zero entries, this col is entirely null, still free variable
        else:
            null_cols.add(j)

    # finds all combinations of null free variables, adds combinations to each basis vector
    all_null = []
    for i in range(1, len(null_cols) + 1):
        all_null += list(combinations(null_cols, i))

    for v in kernel_basis:
        for tup in all_null:
            for c in tup:
                new_basis = kernel_basis[v][0][:]
                new_basis[c] = 1
                kernel_basis[v].append(new_basis)

    return kernel_basis


@deprecated
def _factor_with_known(p, q, n):
    """Helper function for all integer factoring functions, which further factors integer given known factors.

    :param p: integer that divides n, not necessarily prime
    :param q: same as p
    :param n: integer to be factored
    :return: dictionary. keys: all primes factors of n, values: powers of prime factors"""

    # if inputs are non-trivial, try to factor more, if trivial, ignore
    fact_p = quadratic_sieve(p) if p not in [1, n] else {}
    fact_q = quadratic_sieve(q) if q not in [1, n] else {}

    # if failed to factor p or q further, add them to dictionary as non-prime factors
    if fact_q is None:
        fact_q = {q: 1}
    if fact_p is None:
        fact_p = {p: 1}

    factors_known = join_dict(fact_q, fact_p)

    factors = {}
    for f in factors_known:
        factors[f] = 0

    for f in factors:
        while n % f == 0:
            n //= f
            factors[f] += 1

    if n == 1:
        return factors
    if is_prime(n):
        return join_dict(factors, {n: 1})

    more_factors = quadratic_sieve(n)
    return join_dict(factors, {n: 1}) if more_factors is None else join_dict(more_factors, factors)


@deprecated
def __quadratic_sieve1(n, B=None):
    """Third attempt at quadratic sieve. Keeping around for debugging and learning purposes.
    Don't use, wont work!"""

    from math import e

    if B is None:
        L = pow(e, sqrt(log(n) * log(log(n))))
        B = int(pow(L, 1 / sqrt(2))) + 10
    print(f"B: {B}")
    primes = primes_lt(B)

    bases, squares, exp = find_perfect_squares(n, primes)

    # print(exp)
    # print("matrix: ")

    matrix = None

    # print(matrix)
    # print("break")
    m = gaussian_elimination_mod(matrix)
    basis = kernel(m)

    print([basis[7][0]] * matrix)

    def write():
        with open("test.txt", "w") as f:
            string = "{\n"
            for v in basis:
                string += " " + str(v) + ":"
                for row in basis[v]:
                    string += "\n   " + str(row)
                string += "\n\n"
            string += "}"
            f.write(string)

    return None
    width = len(exp[0])
    null_array = numpy.array([0] * width, dtype=object)
    for i in range(2, length := len(exp)):
        indices = [0] * i
        choice = null_array.copy()

        # adds each row of choice together mod 2, if zero vector, this solution possible
        for row in gen_choice(indices, matrix):
            choice = (choice + row) % 2

        # finds all positions in row where value is one, returning positions in tuple containing list
        ones = numpy.where(choice)

        # if there are no instances of row having 1, this row vector represents a perfect square
        if len(ones[0]) == 0:
            a = 1

            # if don't throw dtype=object numpy complains about adding object type to array of 0's
            e = numpy.array([0] * width, dtype=object)
            for j in indices:
                a *= bases[j]
                e += matrix[j]
            # print(e // 2, primes, exp_value(e//2, primes=primes))
            b = exp_value(e // 2, primes=primes)
            if pow(isqrt(b), 2) != b:
                print(a, b, pow(b, 2))
                print(f"exponents: {e}")
                print(f"a: {a}, a^2 - n = {pow(a, 2) - n}")
                print(f"factors of result: {factor_if_smooth(pow(a, 2) - n, primes)}")
                return None
            p, q = gcd(a + b, n), gcd(a - b, n)
            if 1 < p < n or 1 < q < n:
                return _factor_with_known(p, q, n)
            else:

                # debugging step check if taking gcd correctly, if p, q wasn't between 1, n, one of them should be n
                assert n in [p, q]

        indices[-1] += 1

        # calculates the index of rows to be chosen s.t. all possible combinations w/ replacement will be searched
        for index in range(len(indices) - 1, 0, -1):
            j = indices[index]

            # if reached max index or at last index and incrementing would reach max index, reset
            if j == length:
                indices[index - 1] += 1
                indices[index] = 0

        if indices[0] == length:
            return None


@deprecated
def __kernel_old(matrix):
    # dictionary with keys as index of column of each pivot, values as index of each non-zero entry in same
    # row as pivot (ex. if pivot in col 2 of row 1, kernel = {2: [4, 5]}, this is stating that columns
    # 4 and 5 in thw row of this pivot were non-zero and are therefore free variables
    kernel = {}

    # dictionary with keys as row indices of each pivot and values as column indices of each pivot, this is a
    # helper dictionary to assist with construction of basis, has no other use, try to remove this if possible
    pivots = {}

    # set of column indices of all free variables
    free_vars = set()

    # set of column indices of all columns with entirely zero entries
    null_vars = set({i for i in range(len(matrix[0]))})

    # set of column indices of all pivot variables
    pivot_set = set()

    # iterate across matrix row-wise
    for j in range(len(matrix[0])):

        # iterate down columns
        for i in range(len(matrix)):
            e = matrix[i][j]

            # if we hit a non-zero entry and are in a pivot row, add this entry as a free variable to kernel
            # also add column index as index of free variable to set free_vars
            if e == 1 and i in pivots:
                free_vars.add(j)
                kernel[pivots[i]].append(j)

            # otherwise, a non-zero entry is indicative of a pivot, so add to pivot set, pivot dict, and initialize
            # kernel entry to be empty list, which will be filled with indices of free variables
            elif e == 1:
                pivot_set.add(j)
                pivots[i] = j
                kernel[j] = []

    # null_vars initialized to be indices of all columns, now removing indices of all non-zero free variables and
    # indices of all pivots lets null_vars be set of columns with all-zero entries
    null_vars -= free_vars
    null_vars -= pivot_set

    # constructs dictionary with keys as the index of each of the free variables, and values as each of the pivots
    # that are dependent on the specific free variable (ex. {
    kernel_vars = {}

    # iterate through pivots
    for row in pivots:
        col = pivots[row]

        # iterate through list of free variables of each pivot column
        for v in kernel[col]:

            # add to dictionary list of column indices of pivots that depend on the free variable v
            if v in kernel_vars:
                kernel_vars[v].append(col)
            else:
                kernel_vars[v] = [col]

    # constructs basis vectors of kernel, not including null free variables (columns that were entirely null)
    kernel_basis = {}
    for v in kernel_vars:

        # initialize column to be entirely null
        column = [0] * len(matrix[0])
        kernel_basis[v] = [column]

        # make sure that each basis includes reference column (ex. basis [1, 0, 1] * x2 has to also include second row
        # in basis so it should actually be [1, 1, 1])
        kernel_basis[v][0][v] = 1

        # adds list of column indices of pivots that depend on this free variable
        for c in kernel_vars[v]:
            kernel_basis[v][0][c] = 1

    # all_represents a list of all combinations of fully null columns, which are free variables that do not change
    # the list mod 2 (ex. if columns 2 and 5 were null, kernel basis should have its normal vector, a vector where
    # the entry for column 2 is 1, where column 5 is 1, and where both are 1, since neither change the value of
    # the vector mod 2, but when performing calculations in quadratic sieve function, adding these columns does
    # change the result)
    all_null = []
    for i in range(1, len(null_vars) + 1):
        all_null += list(combinations(null_vars, i))

    # iterate through each free var, adding new instances of its basis vector to the list, each a slightly different
    # iteration based off of which null column is present
    for v in kernel_basis:
        for tup in all_null:
            for c in tup:
                new_basis = kernel_basis[v][0][:]
                new_basis[c] = 1
                kernel_basis[v].append(new_basis)

    return kernel_basis
