import os

from math import sqrt, log, gcd
from random import randrange

from cryptography318.factor import lenstra_ecm
from cryptography318.factor import siqs
from cryptography318.prime.prime import isprime, primesieve, next_prime


def factor(n, rho=True, ecm=True, p1=True, qs=True, limit=None):
    """
    Attempts to factor given integer with four methods, returning None if un-factorable.
    Function first checks if number is prime then finds all small factors if any exist.
    Then (assuming no specific methods were set to False) Pollard's Rho, Lenstra's ECM,
    Pollard's P-1, and the Quadratic Sieve are used sequentially to find non-trivial factor
    of n. If any of these methods succeed, a recursive call of factor is used to factor
    the remaining n. All boolean values for use of methods are preserved in the
    recursive call.

    :param n: int number to be factored
    :param rho: bool determining if Pollard's Rho algorithm should be used
    :param ecm: bool determining if Lenstra's ECM algorithm should be used
    :param p1: bool determining if Pollard's P-1 algorithm should be used
    :param qs: bool determining if Quadratic Sieve algorithm should be used
    :param limit: integer limit of factors to be found using small_factors()
    :return: dictionary of all primes factors and their powers, or None if not factorable
    """

    if n == 1:
        return {}
    elif isprime(n):
        return {n: 1}

    if limit is None:
        limit = 32768

    factors = {}
    k, p = factor_small(factors, n, limit)
    if k == 1:
        return factors
    elif k != n:
        if isprime(k):
            factors[k] = factors.get(k, 0) + 1
            return factors
        n = k

    factor_kwargs = {"rho": rho, "ecm": ecm, "p1": p1, "qs": qs, "limit": limit}

    if rho:
        if not _factor_further(n, pollard_rho_factor(n), factors, **factor_kwargs):
            return factors

    if ecm:
        if not _factor_further(n, lenstra_ecm(n), factors, **factor_kwargs):
            return factors

    if p1:
        if not _factor_further(n, pollard_p1(n), factors, **factor_kwargs):
            return factors

    if qs:

        fp = "primes.txt" if os.path.exists("../prime/primes.txt") else None

        # Nothing left after quadratic sieve, so just return factors
        _factor_further(n, siqs(n, fp=fp, loud=False), factors, **factor_kwargs)
        return factors


def pollard_p1(n, B=None, _retry=5):
    """
    Pollard's p - 1 algorithm for factoring large composites.
    Returns one non-trivial factor if factor-able, False if otherwise.

    Pollard's p - 1 is best used to remove smaller factors from a larger composite.
    """

    from math import e

    if B is None:
        L = pow(e, sqrt(log(n) * log(log(n))))
        B = int(pow(L, 1 / sqrt(2)))

    primesieve.extend(B)

    if isprime(n):
        return n

    a = 2
    primes = primesieve[:B]
    for _ in range(_retry):
        m = a
        for j in primes:
            exp = int(log(B, j))
            m = pow(m, pow(j, exp), n)
        q = gcd(m - 1, n)
        if 1 < q < n:
            return q

        a = randrange(2, n - 2)

    return None


def pollard_rho_factor(n, mix=None, _retry=5):
    if n < 10:
        return factor_small({}, n, 10)
    elif not callable(mix):
        def mix(e): return (pow(e, 2, n) + 1) % n

    y = 2
    for _ in range(_retry):
        x = y
        while True:
            x = mix(x)
            y = mix(mix(y))

            q = gcd(abs(x - y), n)
            if q == n:
                break
            elif 1 < q:
                return q

        # If didn't find any, try new mixing function and starting value
        y = randrange(0, n - 1)
        a = randrange(1, n - 3)

        def mix(e): return (pow(e, 2, n) + a) % n

    return None


def factor_small(factors, n, limit):
    """
    Computes all prime factors, up to integer limit, of n when given n is small. Returns
    list of found factors and next odd integer to be checked as factor.
    """

    # Remove as many powers of 2 as possible
    t = 0
    while not n & 1:
        n >>= 1
        t += 1

    if t:
        factors[2] = t

    r = 0
    while (d := divmod(n, 3))[1] == 0:
        n = d[0]
        r += 1

    if r:
        factors[3] = r

    # similarly reduce powers of 'primes' ascending until limit
    p = 5
    while 1:
        r = 0
        while (d := divmod(n, p))[1] == 0:
            n = d[0]
            r += 1

        if r:
            factors[p] = r

        p = next_prime(p)

        # Since all smaller factors have been removed, p is a factor iff kp | n w/ k >= p
        if p > limit or p * p > n:
            break

    return n, p


def _factor_further(n, f, factors, **kwargs):
    """
    Helper function for factor that tries to finish factoring n, having found one
    non-trivial factor f. Function also checks if f needs to be factored. Returns
    False if n has been completely factored, True otherwise.
    """
    if f:
        n //= f
        if isprime(f):
            factors[f] = factors.get(f, 0) + 1
        else:
            factors_f = factor(f, **kwargs)
            for prime in factors_f:
                factors[prime] = factors.get(prime, 0) + factors_f[prime]

        if n == 1:
            return False
        else:
            more_facs = factor(n, **kwargs)
            for prime in more_facs:
                factors[prime] = factors.get(prime, 0) + more_facs[prime]
            return False
    else:
        return True
