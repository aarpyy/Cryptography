from math import isqrt, sqrt, log, gcd
from random import randrange
from sympy.ntheory.primetest import is_square
from functools import reduce

from .prime import isprime, primes_lt_gen, next_prime
from .tools import join_dict
from .elliptic import ecm_mont
from .quadratic_sieve import quadratic_sieve


def factor(n):
    """
    Attempts to factor given integer with four methods, returning None if un-factorable.
    Function first checks if number is prime, then iterates through all primes < 1000 attempting
    to divide n. Function then Lenstra's Elliptic Curve factorization algorithm to try finding small
    factors, then tries Pollard's P-1 algorithm to find one non-trivial factor of n.
    If it succeeds, adds factor to solution set. If n is still factorable, tries quadratic sieve method
    to return all remaining factors. If this returns None, uses sympy's factorint() method and returns result.

    :param n: int number to be factored
    :return: dictionary of all primes factors and their powers, or None if not factorable

    Note
    ----
    It is technically possible for the factors returned to be non-prime factors, but this is only possible
    in situations where n is factorable but the non-prime factors of n are not factorable. Since this is a very
    rare situation, it can be assumed that all returned factors are prime.
    """

    if n == 1:
        return {}

    factors = {}

    n, p = factor_small(factors, n, 1024)

    if n == 1:
        return factors

    if isprime(n):
        return join_dict(factors, {n: 1})

    if is_square(n):
        return _reduce_factors({isqrt(n): 2})

    one_factor = ecm_mont(n)
    n //= one_factor
    factors = join_dict(factors, {one_factor: 1}) if isprime(one_factor) else join_dict(factors, factor(one_factor))
    if isprime(n) or n == 1:
        return factors

    k = pollard_p1(n)
    if isinstance(k, dict):
        return join_dict(factors, k)

    factors_qs = quadratic_sieve(n)
    if factors_qs is None:
        return join_dict(factors, {n: 1})

    factors = join_dict(factors, factors_qs)
    return _reduce_factors(factors)


def pollard_p1(n, B=None, _retry=5):
    """Pollard's p - 1 algorithm for factoring large composites.
    Returns one non-trivial factor if factor-able, False if otherwise."""

    from math import e

    if B is None:
        L = pow(e, sqrt(log(n) * log(log(n))))
        B = int(pow(L, 1 / sqrt(2)))

    if isprime(n):
        return n

    a = 2
    for _ in range(_retry):
        m = a
        for j in primes_lt_gen(B):
            exp = int(log(B, j))
            m = pow(m, pow(j, exp), n)
        q = gcd(m - 1, n)
        if 1 < q < n:
            return _reduce_factors({q: 1, n // q: 1})

        a = randrange(2, n - 2)

    return None


def pollard_rho_factor(n, mix=None, _retry=5):
    if not callable(mix):
        mix = lambda e: (pow(e, 2, n) + 1) % n

    y = 2
    for _ in range(_retry):
        x = y
        while True:
            x = mix(x)
            y = mix(mix(y))

            q = gcd(abs(x - y), n)
            if q == n:
                break
            if 1 < q:
                return _reduce_factors({q: 1, n // q: 1})

        # if didn't find any, try new mixing function and starting value
        y = randrange(0, n - 1)
        a = randrange(1, n - 3)
        mix = lambda e: (pow(e, 2, n) + a) % n

    return None


def _reduce_factors(factors):
    """Given a dictionary of prime or non-prime factors of an integer n, computes all
    prime factors of n and returns in a dictionary, with keys as prime factors
    and values as powers of each prime factor.

    Note
    ----
    Refer to note in factor() about possibility of non-prime factors being returned
    """

    # outer reduce iterates over list of factors of n with an initial value of {}
    return reduce(

        # inner reduce is called iff a factor of n is non-prime; it iterates over factors of non-prime factor of n
        lambda i, c: join_dict(i, reduce(

            # innermost lambda joins dictionary of prime factors of the non-prime factor, multiplying the powers
            # of each prime factor by the power of the original non-prime factor; ex: non-prime factor: {4: 2}
            # would become {2: 2 * {4: 2}[4]} = {2: 2 * 2} = {2: 4}
            lambda a, b: join_dict(a, {b: k[b] * factors[c]}), k, {}

            # if factor of n is prime, join it with the current list of prime factors of n
            # if unable to factor non-prime factor, just add
        )) if (not isprime(c) and (k := factor(c)) is not None) else join_dict(i, {c: factors[c]}), factors, {}
    )


def factor_with_base(n, factors):
    """Function that takes in number and list of factors and returns a list of each of the powers of each factor."""

    exp = [0] * len(factors)
    for i, f in enumerate(factors):
        while n % f == 0:
            n //= f
            exp[i] += 1
    return exp


def trailing(n):
    """
    Computes the number of trailing 'zeros' of given integer in binary
    representation. Equivalent to asking the question: how many times can
    I right shift this integer until it would be converted to a float.
    """
    count = 0
    while not n & 1:
        count += 1
        n >>= 1
    return count


def multiplicity(p, n):
    """
    Computes smallest integer m s.t. p**m | n. Returns 0 if p does not
    divide n.
    """
    b = p
    m = 1
    while not n % b:
        b *= b
        m *= 2

    while n % b:
        b //= p
        m -= 1

    return m


def factor_small(factors, n, limit):
    """
    Computes all prime factors, up to integer limit, of n when given n is small. Returns
    list of found factors and next odd integer to be checked as factor.
    """

    k = n

    # remove as many powers of 2 as possible
    t = trailing(n)
    if t:
        n >>= t
        factors[2] = t

    r = 0
    while not n % 3:
        n //= 3
        r += 1

    if r:
        factors[3] = r

    # similarly reduce powers of 'primes' ascending until limit
    p = 3
    while 1:

        p += 2
        if p > limit:
            break

        r = 0
        while not n % p:
            n //= p
            r += 1

        if r:
            factors[p] = r

        if p * p > k:
            break

        p += 4
        if p > limit:
            break

        r = 0
        while not n % p:
            n //= p
            r += 1

        if r:
            factors[p] = r

        if p * p > k:
            break

    return n, p
