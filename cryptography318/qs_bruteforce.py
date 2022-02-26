from math import sqrt, log, gcd, prod, isqrt
from sympy.ntheory.primetest import is_square
import linalg
from utils import smooth_factor
from prime import isprime, primesieve, randprime


nprimes = 0
primes = []  # primes specifically less than B

required_relations_ratio = 1.05


def binary_dot(a, b):
    return sum(y for x, y in zip(a, b) if x)


def exp_value(exp):
    """Calculates the value of a list of powers of primes. If only p is given, assumes list of primes to be
    from 2 to largest prime <= p. If list of exponents does not match the powers of the continuous ascending
    list of primes, this will compute incorrectly."""

    global primes

    # raises each prime to the corresponding power in list exp, then reduces that list with multiplication
    return prod(map(lambda p, e: pow(p, e), primes, exp))


def find_perfect_squares(n, required_relations):
    """Helper function for Quadratic Sieve that generates N = len(primes) integers that minus p are perfect squares
    and are also B-smooth.

    :param n: number to be factored
    :param required_relations: number of relations needed ot build matrix
    :return: tuple of perfect squares"""

    print(f"required relations: {required_relations}")

    global nprimes

    a = isqrt(n) + 1

    relations_found = 0
    perfect_sq_base = []
    perfect_sq_exp = []

    while relations_found < required_relations:
        b = pow(a, 2) - n
        factors = smooth_factor(b, primes)
        if factors is not None:
            perfect_sq_base.append(a)
            perfect_sq_exp.append(factors)
            relations_found += 1
        a += 1

    return perfect_sq_base, perfect_sq_exp


def quadratic_sieve(n, B=None):
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

    global nprimes, primes
    primesieve.extend(B)
    primes = primesieve[:B]
    nprimes = len(primes)

    print(f"factor base size: {nprimes}")

    # Bases list of a s.t. a^2 is b smooth, exp is powers of b smooth num
    bases, exp = find_perfect_squares(n, int(nprimes * required_relations_ratio))

    mod2 = []
    for i in range(len(exp[0])):
        mod2.append([])
        for j in range(len(exp)):
            mod2[i].append(exp[j][i] % 2)

    kernel = linalg.binary_kernel(mod2)

    for vector in kernel:   # iterate over basis of kernel
        e = map(lambda x: x // 2, linalg.vecmatmul(vector, exp))
        a = 1
        for j, k in zip(vector, bases):
            if j:
                a *= k

        b = exp_value(e)
        p, q = gcd(a + b, n), gcd(a - b, n)

        if 1 < p < n:
            return p
        if 1 < q < n:
            return q

    # this return statement should never hit, if it does consider adding more rows to matrix in function find_perf_sq
    return None


if __name__ == "__main__":
    lower = pow(10, 6)
    upper = pow(10, 7)
    A = randprime(lower, upper)
    B = randprime(lower, upper)
    print(f"{A} * {B} = {(N := A * B)}")
    factor = quadratic_sieve(N, 1000)
    if factor is not None:
        print(f"factor of n: {factor}; n / factor = {N // factor}")
