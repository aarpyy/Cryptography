from math import sqrt, log, gcd, prod, isqrt
from sympy.ntheory.primetest import is_square
from .prime import isprime, primes_lt
from .linalg import dot, matrix, Array, BinaryMatrix


search_time = None


def exp_value(exp, primes):
    """Calculates the value of a list of powers of primes. If only p is given, assumes list of primes to be
    from 2 to largest prime <= p. If list of exponents does not match the powers of the continuous ascending
    list of primes, this will compute incorrectly."""

    # raises each prime to the corresponding power in list exp, then reduces that list with multiplication
    return prod(map(lambda p, e: pow(p, e), primes, exp))


def factor_if_smooth(n, primes):
    """Helper function for quadratic sieve that returns factors of n if n can be factored
    only by list of given primes, returning None otherwise."""

    exp = [0] * len(primes)

    for i, p in enumerate(primes):
        while n % p == 0:
            n //= p
            exp[i] += 1

    return exp if n == 1 else None


def find_perfect_squares(n, primes):
    """Helper function for Quadratic Sieve that generates N = len(primes) integers that minus p are perfect squares
    and are also B-smooth.

    :param n: int
    :param primes: list of primes
    :return: tuple of perfect squares"""

    # attempt to not have to generate a lot of extra rows for smaller n, tweak this
    extra_rows = len(primes) // 10

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
        while i < len(primes) + extra_rows:
            b = pow(a, 2) - n
            factors = factor_if_smooth(b, primes)
            if factors is not None:
                perfect_sq_base.append(a)
                perfect_sq_exp.append(Array(factors))
                i += 1
                curr = time()

            # if too hard to find b smooth nums just try with already found, needs to be at least square matrix
            if i >= len(primes) and time() - curr > search_time:
                return perfect_sq_base, perfect_sq_exp
            a += 1
    else:
        while i < len(primes) + extra_rows:
            b = pow(a, 2) - n
            factors = factor_if_smooth(b, primes)
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

    primes = primes_lt(B)

    global search_time
    search_time = force

    bases, exp = find_perfect_squares(n, primes)  # bases list of a s.t. a^2 is b smooth, exp is powers of b smooth num

    mat = matrix(exp).transpose()  # mat needs to be transpose and kernel also wants transpose so do this step now

    basis = BinaryMatrix(
        list(map(lambda r: r.mod(2), mat.array))  # use BinaryMatrix instead of matrix to reduce parse time of object
    ).kernel().transpose()  # basis of kernel as matrix

    for v in basis:  # iterate over basis of kernel
        e = map(lambda x: x // 2, map(lambda y: dot(v, y), mat))
        a = dot(v, bases)

        b = exp_value(e, primes)
        p, q = gcd(a + b, n), gcd(a - b, n)

        # if one solution is non-trivial factor, _factor_with_known will take care of the rest
        if 1 < p < n:
            return {p: 1, n//p: 1}
        if 1 < q < n:
            return {q: 1, n//q: 1}

    # this return statement should never hit, if it does consider adding more rows to matrix in function find_perf_sq
    return None
