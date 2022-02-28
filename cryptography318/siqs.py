import random
import sys

from math import log2, sqrt, isqrt, gcd, ceil
from functools import reduce

from utils import n_digits, smooth_factor, exp_value
from linalg import binary_kernel, vecmatmul
from prime import primesieve, quadratic_residue, sqrt_mod, randprime


min_factor = 2000
max_factor = 4000
min_n_factors = 20
trials_a = 30
trial_error_margin = 25
required_relations_ratio = 1.05

relations_found = 0
min_sieve = 0

a = b = 0               # type: int
B = []                  # type: list[int]
primes = []             # type: list[int]
t_sqrt = []             # type: list[int]
log_p = []              # type: list[int]
factor_base = []        # type: list[int]
soln1 = []              # type: list[int]
soln2 = []              # type: list[int]
B_ainv_2 = []           # type: list[list[int]]
a_factors = set()       # type: set[int]
a_non_factors = set()   # type: set[int]
smooth_u = []           # type: list[list[int]]
smooth_t = []           # type: list[int]

rand = random.Random()


class QSPoly:

    __slots__ = "args",

    def __new__(cls, *args):
        self = super(QSPoly, cls).__new__(cls)
        self.args = [*args]
        return self

    def __call__(self, x: int) -> int:
        return reduce(lambda y, z: (y * x) + z, self.args, 0)


def choose_f(digits: int):
    if digits < 38:
        return 4200
    elif digits < 40:
        return 5600
    elif digits < 42:
        return 7000
    elif digits < 44:
        return 8400
    elif digits < 48:
        return 13000
    elif digits < 52:
        return 16000
    elif digits < 56:
        return 29000
    elif digits < 60:
        return 60000
    elif digits < 70:
        return 100000
    elif digits < 80:
        return 350000
    else:
        return 900000


def choose_m(digits: int):
    if digits < 45:
        return 40000
    elif digits < 52:
        return 65536
    elif digits < 88:
        return 196608
    else:
        return 589824


def init_siqs(n, *, file=None):

    global factor_base, t_sqrt, log_p, primes, soln1, soln2

    F = choose_f(n_digits(n))

    print(f"F: {F}")

    p = 1
    if file is None:
        primesieve.extend(F)
        primes = primesieve[:F]
    else:
        with open(file, "r") as prime_file:
            while p < F:
                primes.append(p := int(prime_file.readline()))

    print(f"primes < F: {len(primes)}")

    for prime in primes:
        if quadratic_residue(n, prime):
            factor_base.append(prime)

            t_sqrt.append(sqrt_mod(n, prime))
            log_p.append(round(log2(prime)))

    print(f"Size of factor base: {len(factor_base)}")

    soln1 = [None] * len(factor_base)
    soln2 = [None] * len(factor_base)


def smooth_a(n, m):
    global a, factor_base, a_factors, a_non_factors

    s = len(factor_base)

    start = 0

    while factor_base[start] < min_factor:
        start += 1
        if start >= s:
            start = 0
            break

    stop = start
    while factor_base[stop] < max_factor:
        stop += 1
        if stop >= s:
            stop = s - 1
            break

    if stop - start < min_n_factors:
        raise ValueError("Not enough factors in factor base, try increasing F")

    target = isqrt(n + n) // m
    min_a = target / sqrt((factor_base[stop] + factor_base[start]) / 2)
    opt_ratio = 0.9

    a_factors = set()
    best_ratio = None   # type: None | float

    for _ in range(trials_a):

        A = 1
        tmp_factors = set()
        while A < min_a:
            i = rand.randrange(start, stop)
            if i not in tmp_factors:
                tmp_factors.add(i)
                A *= factor_base[i]

        ratio = A / target
        if best_ratio is None or best_ratio > ratio >= opt_ratio or ratio >= opt_ratio > best_ratio:
            best_ratio = ratio
            a = A
            a_factors = tmp_factors

    set_fb = set(range(len(factor_base)))
    a_non_factors = set_fb - a_factors


def first_poly(n, m):

    global a, b, B, factor_base, t_sqrt, B_ainv_2, a_factors, a_non_factors

    smooth_a(n, m)

    B = []
    sorted_factors = sorted(list(a_factors))
    n_factors = len(a_factors)

    for j in range(n_factors):
        q = factor_base[sorted_factors[j]]

        assert a % q == 0

        a_l = a // q
        gamma = (t_sqrt[j] * pow(a_l, -1, q)) % q
        if gamma > q / 2:
            gamma = q - gamma

        B.append(a_l * gamma)

    b = sum(B) % a

    _b = a - b if b + b > a else b

    size_fb = len(factor_base)
    B_ainv_2 = []
    for _ in range(len(a_factors)):
        B_ainv_2.append([None] * size_fb)

    for p in a_non_factors:
        prime = factor_base[p]
        a_inv = pow(a, -1, prime)
        for j in range(n_factors):
            B_ainv_2[j][p] = ((B[j] + B[j]) * a_inv) % prime

        t = t_sqrt[p]
        soln1[p] = (a_inv * (t - b)) % prime
        soln2[p] = (a_inv * (-t - b)) % prime

    return QSPoly(a * a, 2 * a * _b, _b * _b - n), QSPoly(a, _b)


def next_poly(i, n):

    global b, B, a_non_factors, soln1, soln2, B_ainv_2, factor_base

    v = 1
    j = i
    while not j & 1:
        j >>= 1
        v += 1

    sign = -1 if ceil(i / pow(2, v)) & 1 else 1

    v -= 1

    b = (b + 2 * sign * B[v]) % a
    _b = a - b if b + b > a else b

    for p in a_non_factors:
        prime = factor_base[p]
        soln1[p] = (soln1[p] + sign * B_ainv_2[v][p]) % prime
        soln2[p] = (soln2[p] + sign * B_ainv_2[v][p]) % prime

    return QSPoly(a * a, 2 * a * _b, _b * _b - n), QSPoly(a, _b)


def sieve(m):
    global a_non_factors, factor_base, log_p, soln1, soln2

    m2_1 = m + m + 1
    sieve_array = [0] * m2_1

    for p in a_non_factors:
        prime = factor_base[p]

        i_min = -((m + soln1[p]) // prime)
        for j in range(soln1[p] + m + i_min * prime, m2_1, prime):
            sieve_array[j] += log_p[p]

        if prime != 2:
            i_min = -((m + soln2[p]) // prime)
            for j in range(soln2[p] + m + i_min * prime, m2_1, prime):
                sieve_array[j] += log_p[p]

    return sieve_array


def trial_division(sieve_array, m, g: QSPoly, h: QSPoly):
    global primes, smooth_t, smooth_u, relations_found, min_sieve

    for i, s in enumerate(sieve_array):
        if s >= min_sieve:
            x = i - m
            u = g(x)
            if (powers := smooth_factor(u, primes)) is not None:
                t = h(x)
                smooth_u.append(powers)
                smooth_t.append(t)
                relations_found += 1


def solve_matrix(n):
    global smooth_t, smooth_u, primes

    mod2 = []
    for i in range(len(smooth_u[0])):
        mod2.append([])
        for j in range(len(smooth_u)):
            mod2[i].append(smooth_u[j][i] % 2)

    kernel = binary_kernel(mod2)

    for vector in kernel:  # iterate over basis of kernel
        e = map(lambda v: v // 2, vecmatmul(vector, smooth_u))
        x = 1
        for j, k in zip(vector, smooth_t):
            if j:
                x *= k

        y = exp_value(e, primes)
        p, q = gcd(x + y, n), gcd(x - y, n)

        if 1 < p < n:
            return p
        if 1 < q < n:
            return q

    return None


def main(n):
    global min_sieve, factor_base, relations_found, a_factors

    init_siqs(n, file="primes.txt")

    m = choose_m(n_digits(n))
    min_sieve = log2(isqrt(n) * m) - trial_error_margin

    g, h = first_poly(n, m)

    required_relations = int(len(factor_base) * required_relations_ratio)

    n_poly = pow(2, len(a_factors) - 1) - 1
    i = 1
    while relations_found < required_relations:
        sieve_array = sieve(m)
        trial_division(sieve_array, m, g, h)

        if i >= n_poly:
            g, h = first_poly(n, m)
            n_poly = pow(2, len(a_factors) - 1) - 1
            i = 0
            # print(f"Relations found: {relations_found}/{required_relations}")
        else:
            g, h = next_poly(i, n)

        i += 1

    if (factor := solve_matrix(n)) is not None:
        print(f"Factor: {factor}")


if __name__ == "__main__":
    lower = pow(10, 15)
    upper = pow(10, 16)
    X = randprime(lower, upper)
    Y = randprime(lower, upper)
    print(f"{X} * {Y} = {(N := X * Y)}")
    main(N)

