import random
import sys

from math import log2, sqrt, isqrt, prod

from utils import n_digits
from prime import primesieve, quadratic_residue, quadratic_non_residue, sqrt_mod


min_factor = 2000
max_factor = 4000
min_n_factors = 20
trials_a = 30

a = b = 0

rand = random.Random()


class QSPoly:

    __slots__ = "coeffs",

    def __new__(cls, *args):
        self = super(QSPoly, cls).__new__(cls)
        self.coeffs = [*args]
        return self

    def __call__(self, x: int) -> int:
        acc = 0
        for a in self.coeffs:
            acc = (acc * x) + a
        return acc


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
    F = choose_f(n_digits(n))

    primes = []
    p = 1
    if file is None:
        primesieve.extend(F)
        primes = primesieve[:F]
    else:
        with open(file, "r") as prime_file:
            while p < F:
                primes.append(p := int(prime_file.readline()))

    t_sqrt = []         # type: list[int]
    log_p = []          # type: list[int]
    factor_base = []    # type: list[int]

    for prime in primes:
        if quadratic_residue(n, prime):
            factor_base.append(prime)

            t_sqrt.append(sqrt_mod(n, prime))
            log_p.append(round(log2(prime)))

    return factor_base, t_sqrt, log_p


def smooth_a(n, m, fb):
    global a

    s = len(fb)

    start = 0

    while fb[start] < min_factor:
        start += 1
        if start >= s:
            start = 0
            break

    stop = start
    while fb[stop] < max_factor:
        stop += 1
        if stop >= s:
            stop = s - 1
            break

    if stop - start < min_n_factors:
        raise ValueError("Not enough factors in factor base, try increasing F")

    target = isqrt(n + n) // m
    min_a = target / sqrt((fb[stop] + fb[start]) / 2)
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
                A *= fb[i]

        ratio = A / target
        if best_ratio is None or best_ratio > ratio >= opt_ratio or ratio >= opt_ratio > best_ratio:
            best_ratio = ratio
            a = A
            a_factors = tmp_factors

    set_fb = set(fb)
    a_non_factors = set_fb - a_factors
    return a_factors, a_non_factors


def first_poly(n, m, fb, t_sqrt, soln1, soln2, B_ainv_2):

    global a, b

    a_factors, a_non_factors = smooth_a(n, m, fb)

    B = []
    sorted_factors = sorted(list(a_factors))
    n_factors = len(a_factors)

    for j in range(n_factors):
        q = fb[sorted_factors[j]]

        assert a % q == 0

        a_l = a // q
        gamma = (t_sqrt[j] * pow(a_l, -1, q)) % q
        if gamma > q / 2:
            gamma = q - gamma

        B.append(a_l * gamma)

    b = sum(B) % a

    _b = a - b if b + b > a else b

    for p in a_non_factors:
        prime = fb[p]
        a_inv = pow(a, -1, prime)
        for j in range(n_factors):
            B_ainv_2[j][p] = ((B[j] + B[j]) * a_inv) % prime

        t = t_sqrt[p]
        soln1[p] = (a_inv * (t - b)) % prime
        soln2[p] = (a_inv * (-t - b)) % prime

    return QSPoly(a * a, 2 * a * _b, _b * _b - n), QSPoly(a, _b)


def main(n):
    m = choose_m(n_digits(n))
    factor_base, t_sqrt, log_p = init_siqs(n)
    a_factors, a_non_factors = smooth_a(n, m, factor_base)
    size = len(factor_base)
    B_ainv_2 = []
    for _ in range(len(a_factors)):
        B_ainv_2.append([0] * size)

    soln1 = []
    soln2 = []

    g, h = first_poly(n, m, factor_base, t_sqrt, soln1, soln2, B_ainv_2)


if __name__ == "__main__":
    main(sys.argv[0])

