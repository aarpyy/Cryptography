from math import sqrt, log, gcd, prod, isqrt
from sympy.ntheory.primetest import is_square
from .prime import isprime, multiplicity
from .linalg import dot, matrix, Array, BinaryMatrix


search_time = None
nprimes = 0
primes = []  # primes specifically less than B


def exp_value(exp):
    """Calculates the value of a list of powers of primes. If only p is given, assumes list of primes to be
    from 2 to largest prime <= p. If list of exponents does not match the powers of the continuous ascending
    list of primes, this will compute incorrectly."""

    global primes

    # raises each prime to the corresponding power in list exp, then reduces that list with multiplication
    return prod(map(lambda p, e: pow(p, e), primes, exp))


def factor_if_smooth(n):
    """Helper function for quadratic sieve that returns factors of n if n can be factored
    only by list of given primes, returning None otherwise."""

    global primes, nprimes

    exp = [0] * nprimes

    for i, p in enumerate(primes):
        r = 0
        while n % p == 0:
            n //= p
            r += 1
            if r >= 20:
                rr = multiplicity(p, n)
                n //= p ** rr
                r += rr
        if r:
            exp[i] = r

    return exp if n == 1 else None


def find_perfect_squares(n):
    """Helper function for Quadratic Sieve that generates N = len(primes) integers that minus p are perfect squares
    and are also B-smooth.

    :param n: int
    :return: tuple of perfect squares"""

    global nprimes
    # attempt to not have to generate a lot of extra rows for smaller n, tweak this
    extra_rows = isqrt(nprimes)

    a = isqrt(n) + 1

    i = 0
    perfect_sq_base = []
    perfect_sq_exp = []

    # finds + extra_rows so that resulting matrix will have more rows than columns which increases the chances of each
    # row being linearly dependent mod 2
    global search_time

    if isinstance(search_time, int):
        from time import time
        curr = time()
        while i < nprimes + extra_rows:
            b = pow(a, 2) - n
            factors = factor_if_smooth(b)
            if factors is not None:
                perfect_sq_base.append(a)
                perfect_sq_exp.append(Array(factors))
                i += 1
                curr = time()

            # if too hard to find b smooth nums just try with already found, needs to be at least square matrix
            if time() - curr > search_time:
                return perfect_sq_base, perfect_sq_exp
            a += 1
    else:
        while i < nprimes + extra_rows:
            b = pow(a, 2) - n
            factors = factor_if_smooth(b)
            if factors is not None:
                perfect_sq_base.append(a)
                perfect_sq_exp.append(Array(factors))
                i += 1
            a += 1

    return perfect_sq_base, perfect_sq_exp


def quadratic_sieve(n, B=None, force=60):
    """Performs the single polynomial quadratic sieve.
    Attempts to factor given integer n through finding B-smooth perfect square solutions to the
    quadratic function a^2 - n for √n < a with the goal of finding a solution that allows a multiple
    of a factor n to be found, thus allowing the gcd(solution, n) to yield a non-trivial solution.
    Force parameter set to 60s initially, meaning for each b-smooth number 60s is allowed between
    finding each number, if 60s is reached the values found will be returned regardless if found
    enough values for square matrix.

    :param n: integer to be factored
    :param B: required smoothness of solution (all factors of solution have to be <= B), it is recommended to leave
    this value at its default
    :param force: integer value of seconds allowed for searching for each b-smooth number, default 60s
    :return: dictionary of all prime factors of n and their powers, or None if n not-factorable"""

    from math import e

    if n == 1:
        return {}

    if isprime(n):
        return {n: 1}

    if is_square(n):
        return {isqrt(n): 2}

    if B is None:
        L = pow(e, sqrt(log(n) * log(log(n))))
        B = int(pow(L, 1 / sqrt(2)))

    global nprimes, primesieve, primes
    primesieve.extend(B)
    primes = primesieve[:B]
    nprimes = len(primes)

    global search_time
    search_time = force

    bases, exp = find_perfect_squares(n)  # bases list of a s.t. a^2 is b smooth, exp is powers of b smooth num

    mat = matrix(exp).transpose()  # mat needs to be transpose and kernel also wants transpose so do this step now

    basis = BinaryMatrix(
        list(map(lambda r: r.mod(2), mat.array))  # use BinaryMatrix instead of matrix to reduce parse time of object
    ).kernel().transpose()  # basis of kernel as matrix

    for v in basis:  # iterate over basis of kernel
        e = map(lambda x: x // 2, map(lambda y: dot(v, y), mat))
        a = dot(v, bases)

        b = exp_value(e)
        p, q = gcd(a + b, n), gcd(a - b, n)

        if 1 < p < n:
            return p
        if 1 < q < n:
            return q

    # this return statement should never hit, if it does consider adding more rows to matrix in function find_perf_sq
    return None
