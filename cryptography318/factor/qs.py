from collections.abc import Callable
from dataclasses import dataclass
from functools import reduce
from random import Random

import numpy as np
from math import ceil, gcd, isqrt, log2, prod, sqrt

from cryptography318.linalg import kernel_gf2
from cryptography318.prime.prime import quadratic_residue, sqrt_mod, primesieve
from utils.siqs import choose_f, n_digits, choose_m


@dataclass
class FactorBaseItem:
    prime: int
    soln1: int | None
    soln2: int | None
    log_p: int
    B_ainv_2: list[int] | None
    t_sqrt: int
    a_inv: int | None


# Constants for SIQS algorithm
MIN_A_FACTOR = 2000  # Smallest factor of a
MAX_A_FACTOR = 4000  # Largest factor of a
MIN_N_FACTORS = 20
TRIALS_A = 30
TRIAL_ERROR_MARGIN = 25
REQUIRED_RELATIONS_RATIO = 1.05  # Required relations-found:factor-base-size ratio (how tall should matrix be)
TRIALS_LINALG = 5  # Number of allowed attempts at solving linear system before giving up

rand = Random()

class QSPoly(Callable[[int], int]):
    __slots__ = "args", "a", "b"

    def __new__(cls, *args, a=None, b=None):
        self = super().__new__(cls)
        self.args = [*args]
        self.a = a
        self.b = b
        return self

    def __call__(self, x):
        return reduce(lambda y, z: (y * x) + z, self.args, 0)


def vec_matmul_T(vector, matrix):
    """
    Vector x matrix multiplication that takes in matrix
    already transposed, to save time when multiplying against the
    same matrix repeatedly.
    """

    # Cast to Python int so that when we take power later it doesn't throw overflow error
    return (int(vector @ row) // 2 for row in matrix)


class SIQS(object):

    def __init__(self, n, F, *, fp=None, verbose=False):
        self.verbose = verbose
        self.n = n
        self.smooth_t = []
        self.smooth_u = []

        # Choose F based on size of n base 10
        F = F or choose_f(n_digits(n))
        self.print(f"F: {F}")

        # Load primes from file is one was provided
        if fp is not None:
            primesieve.load(fp)

        # Initialize factor base to primes for which n is a quadratic residue
        self.factor_base = []
        for prime in primesieve.primerange(F):
            if quadratic_residue(n, prime):
                self.factor_base.append(FactorBaseItem(prime, None, None, round(log2(prime)), None,
                                                       sqrt_mod(n, prime), None))

        self.print(f"Size of factor base: {len(self.factor_base)}")

        self.required_relations = int(len(self.factor_base) * REQUIRED_RELATIONS_RATIO)

        # Choose sieve range and minimum sieve value
        self.m = choose_m(n_digits(n))

    def print(self, *args, **kwargs):
        if self.verbose:
            print(*args, **kwargs)

    def smooth_a(self):
        """
        Computes and returns coefficient a that is the product of several primes, ideally
        between 2000 and 4000, all in the factor base of primes that n has a quadratic residue
        for.
        """

        s = len(self.factor_base)

        start = 0

        while self.factor_base[start].prime < MIN_A_FACTOR:
            start += 1
            if start >= s:
                start = 0
                break

        stop = start
        while self.factor_base[stop].prime < MAX_A_FACTOR:
            stop += 1
            if stop >= s:
                stop = s - 1
                break

        if stop - start < MIN_N_FACTORS:
            raise ValueError("Not enough factors in factor base, try increasing F")

        target = isqrt(self.n + self.n) // self.m
        min_a = target / sqrt((self.factor_base[stop].prime + self.factor_base[start].prime) / 2)
        opt_ratio = 0.9

        best_ratio = None  # type: None | float
        best_a = None
        best_factors = None
        # Try several ones to find the approximately closest to our target
        for _ in range(TRIALS_A):

            A = 1
            tmp_factors = set()
            while A < min_a:
                i = rand.randrange(start, stop)
                if i not in tmp_factors:
                    tmp_factors.add(i)
                    A *= self.factor_base[i].prime

            ratio = A / target
            if best_ratio is None or best_ratio > ratio >= opt_ratio or ratio >= opt_ratio > best_ratio:
                best_ratio = ratio
                best_a = A
                best_factors = tmp_factors

        return best_a, best_factors

    def first_poly(self):
        """
        Given number to be factored and sieve range, compute `a` as the product of primes in the
        factor base, and from that b such that a | b * b - n. Use these coefficients to
        create two polynomials, one used for finding smooth numbers and the other for finding
        the square root of the value square to find a smooth output.
        :return: two polynomials (ax + b)^2 - n and ax + b
        """

        a, factors = self.smooth_a()

        B = []

        for q in sorted(factors):
            factor = self.factor_base[q]
            a_l = a // factor.prime
            gamma = (factor.t_sqrt * pow(a_l, -1, factor.prime)) % factor.prime
            if gamma > factor.prime / 2:
                gamma = factor.prime - gamma

            B.append(a_l * gamma)

        b = sum(B) % a

        _b = a - b if b + b > a else b

        # assert (_b * _b - self.n) % a == 0

        for factor in self.factor_base:
            # We only want non-factors
            if a % factor.prime == 0:
                continue

            factor.a_inv = pow(a, -1, factor.prime)
            factor.B_ainv_2 = [2 * j * factor.a_inv % factor.prime for j in B]
            factor.soln1 = (factor.a_inv * (factor.t_sqrt - b)) % factor.prime
            factor.soln2 = (factor.a_inv * (-factor.t_sqrt - b)) % factor.prime

        return QSPoly(a * a, 2 * a * _b, _b * _b - self.n, a=a, b=_b), B

    def next_poly(self, i, g, B):
        v = 1
        j = i
        while not j & 1:
            j >>= 1
            v += 1

        sign = -1 if ceil(i / pow(2, v)) & 1 else 1
        v -= 1

        a, b = g.a, g.b
        b = (b + 2 * sign * B[v]) % a
        _b = a - b if b + b > a else b

        # assert (_b * _b - self.n) % a == 0

        for factor in self.factor_base:
            # We only want non-factors
            if a % factor.prime == 0:
                continue

            factor.soln1 = (factor.soln1 + sign * factor.B_ainv_2[v]) % factor.prime
            factor.soln2 = (factor.soln2 + sign * factor.B_ainv_2[v]) % factor.prime

        return QSPoly(a * a, 2 * a * _b, _b * _b - self.n, a=a, b=_b)

    def sieve(self):
        m2 = self.m + self.m
        sieve_array = [0] * (m2 + 1)

        for factor in self.factor_base:
            # We only want non-factors
            if factor.soln1 is None:
                continue

            for j in range((factor.soln1 + self.m) % factor.prime, m2, factor.prime):
                sieve_array[j] += factor.log_p

            if factor.prime == 2:
                continue

            for j in range((factor.soln2 + self.m) % factor.prime, m2, factor.prime):
                sieve_array[j] += factor.log_p

        return sieve_array

    def smooth_factor(self, x):
        exp = [0] * len(self.factor_base)
        for i, f in enumerate(self.factor_base):
            while (d := divmod(x, f.prime))[1] == 0:
                x = d[0]
                exp[i] += 1

        if abs(x) == 1:
            return exp

        return None

    def trial_division(self, sieve_array, g):
        relations_found = 0
        min_sieve = log2(isqrt(self.n) * self.m) - TRIAL_ERROR_MARGIN
        for i, s in enumerate(sieve_array):
            if s < min_sieve:
                continue

            x = i - self.m
            u = g(x)
            if (powers := self.smooth_factor(u)) is not None:
                t = g.a * x + g.b
                self.smooth_u.append(powers)
                self.smooth_t.append(t)
                relations_found += 1

        return relations_found

    def solve_matrix(self):
        T = np.array(self.smooth_u, dtype=object).transpose()
        kernel = kernel_gf2(T)

        for vector in kernel:
            powers = vec_matmul_T(vector, T)

            # We need to do this instead of dot because dot will return the sum instead of product
            x = prod(k for j, k in zip(vector, self.smooth_t) if j)
            y = prod(pow(p.prime, e) for p, e in zip(self.factor_base, powers))
            p, q = gcd(x + y, self.n), gcd(x - y, self.n)

            if 1 < p < self.n:
                return p
            if 1 < q < self.n:
                return q

        return None

    def clear(self):
        self.smooth_u.clear()
        self.smooth_t.clear()


def qs(n, F=None, *, fp=None, verbose=True):
    """
    Performs the Self-Initializing Quadratic Sieve on integer n. For detailed explanation
    of algorithm and sources refer to the full project done in Java at https://github.com/aarpyy/SIQS.
    All references used are linked here but specifics of where these references are used
    can be found in the Java version, as this project is directly adapted from that.

    References
    ----------
    https://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=53C827A542A8A950780D34E79261FF99?doi=10.1.1.26.6924&rep=rep1&type=pdf
    https://github.com/skollmann/PyFactorise/blob/master/factorise.py
    https://www.rieselprime.de/ziki/Self-initializing_quadratic_sieve

    :param n: number to be factored
    :param F: upper bound for primes to be used in factor base
    :param fp: file containing list of primes
    :param verbose: if information should be printed during execution
    :return: factor of n if one exists, otherwise None
    """

    # Initialize factor base, square root N mod p, and log p for all primes p
    # where N is a quadratic residue mod p
    siqs = SIQS(n, F, fp=fp, verbose=verbose)

    # Number of polynomials that can be used with this 'a' value
    i = 1
    last_printed = 0
    relations_found = 0

    # Initialize first polynomial
    g, B = siqs.first_poly()
    for trial in range(TRIALS_LINALG):
        siqs.print("Finding relations...")
        while relations_found < siqs.required_relations:
            sieve_array = siqs.sieve()
            relations_found += siqs.trial_division(sieve_array, g)

            # If we have found anymore relations, print
            if relations_found >= last_printed:
                last_printed = relations_found
                siqs.print(f"\r{relations_found}/{siqs.required_relations}", end="")

            if i >= pow(2, len(B) - 1):
                g, B = siqs.first_poly()
                i = 0
            else:
                g = siqs.next_poly(i, g, B)

            i += 1

        siqs.print()

        if (factor := siqs.solve_matrix()) is not None:
            return factor
        else:
            # Reset relations found and smooth values found
            last_printed = relations_found = 0
            siqs.clear()

    return None
