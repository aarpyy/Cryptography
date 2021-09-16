from math import isqrt, sqrt, log, gcd
from random import randrange
from sympy.ntheory.primetest import is_square

from cryptography318.prime import isprime, primesieve, multiplicity, trailing
from .elliptic import lenstra_ecm
from cryptography318.numbers.quadratic_sieve import quadratic_sieve


def factor(n, rho=True, ecm=True, p1=True, qs=True, limit=None, qs_limit=None):
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
    :param qs_limit: integer time limit of search time for b-smooth numbers in qs algorithm
    :return: dictionary of all primes factors and their powers, or None if not factorable

    Note
    ----
    It is technically possible for the factors returned to be non-prime factors, but this is only possible
    in situations where n is factorable but the non-prime factors of n are not factorable. Since this is a very
    rare situation, it can be assumed that all returned factors are prime.
    """

    if n == 1:
        return {}

    if isprime(n):
        return {n: 1}

    factors = {}
    if limit is None:
        limit = 32768

    k, p = factor_small(factors, n, limit)

    # if factor_small didn't find anything, try to check for square before moving to other algorithms
    if k == n:
        if is_square(n):
            return {isqrt(n): 2}
    n = k

    try:
        _check_factored(factors, n)

        if rho:
            f = pollard_rho_factor(n)
            if f:
                n //= f
                factors[f] = factors.get(f, 0) + 1
                _check_factored(factors, n)  # check to see if more factoring is required

                more_facs = factor(n, rho=rho, ecm=ecm, p1=p1, qs=qs, limit=limit, qs_limit=qs_limit)
                if more_facs:
                    for f, e in more_facs.items():
                        factors[f] = factors.get(f, 0) + e
                    raise StopIteration

        if ecm:
            f = lenstra_ecm(n)
            if f:
                n //= f
                factors[f] = factors.get(f, 0) + 1
                _check_factored(factors, n)  # check to see if more factoring is required

                more_facs = factor(n, rho=rho, ecm=ecm, p1=p1, qs=qs, limit=limit, qs_limit=qs_limit)
                if more_facs:
                    for f, e in more_facs.items():
                        factors[f] = factors.get(f, 0) + e
                    raise StopIteration

        if p1:
            f = pollard_p1(n)
            if f:
                n //= f
                factors[f] = factors.get(f, 0) + 1
                _check_factored(factors, n)  # check to see if more factoring is required

                more_facs = factor(n, rho=rho, ecm=ecm, p1=p1, qs=qs, limit=limit, qs_limit=qs_limit)
                if more_facs:
                    for f, e in more_facs.items():
                        factors[f] = factors.get(f, 0) + e
                    raise StopIteration

        if qs:
            # if any of pollard p-1, pollard rho, or lenstra ecm found a factor then the number is factorable
            # without use of qs since it was able to factor larger N, this should only hit if nothing else can
            # find a factor
            f = quadratic_sieve(n, force=qs_limit)
            if f:
                n //= f
                factors[f] = factors.get(f, 0) + 1
                _check_factored(factors, n)  # check to see if more factoring is required

                more_facs = factor(n, rho=rho, ecm=ecm, p1=p1, qs=qs, limit=limit, qs_limit=qs_limit)
                if more_facs:
                    for f, e in more_facs.items():
                        factors[f] = factors.get(f, 0) + e
                    raise StopIteration

        factors[n] = factors.get(n, 0) + 1  # if this hits then no method could find any factor so just add n to factors

    except StopIteration:
        return _reduce_factors(factors)
    else:
        # if this else hits then no _check_factored raised an error and it can be assumed n was not completely factored
        return factors


def pollard_p1(n, B=None, _retry=5):
    """Pollard's p - 1 algorithm for factoring large composites.
    Returns one non-trivial factor if factor-able, False if otherwise."""

    from math import e

    if B is None:
        L = pow(e, sqrt(log(n) * log(log(n))))
        B = int(pow(L, 1 / sqrt(2)))

    primesieve.extend(B)

    if isprime(n):
        return n

    a = 2
    for _ in range(_retry):
        m = a
        for j in primesieve.range(B):
            exp = int(log(B, j))
            m = pow(m, pow(j, exp), n)
        q = gcd(m - 1, n)
        if 1 < q < n:
            return q

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
                return q

        # if didn't find any, try new mixing function and starting value
        y = randrange(0, n - 1)
        a = randrange(1, n - 3)
        mix = lambda e: (pow(e, 2, n) + a) % n

    return None


def _reduce_factors(factors):
    """Given a dictionary of prime or non-prime factors of an integer n, computes all
    prime factors of n and returns in a dictionary, with keys as prime factors
    and values as powers of each prime factor.

    Notes
    ----
    [1] Refer to note in factor() about possibility of non-prime factors being returned
    [2] Factors is both returned and directly modified, allowing for _reduce_factors() to be
    called without needing to handle return
    """

    new_factors = []  # prime factors of found non-prime factors
    to_del = set()  # set of non-prime factors to be removed
    for f in factors:
        if not isprime(f):
            k = factor(f)
            for a in k:
                k[a] *= factors[f]
            new_factors.append(k)
            to_del.add(f)
    for f in to_del:
        del factors[f]
    for d in new_factors:
        for f in d:
            if f in factors:
                factors[f] += d[f]
            else:
                factors[f] = d[f]
    return factors


def factor_with_base(n, factors):
    """Function that takes in number and list of factors and returns a list of each of the powers of each factor."""

    exp = [0] * len(factors)
    for i, f in enumerate(factors):
        while n % f == 0:
            n //= f
            exp[i] += 1
    return exp


def factor_small(factors, n, limit):
    """
    Computes all prime factors, up to integer limit, of n when given n is small. Returns
    list of found factors and next odd integer to be checked as factor.

    Notes
    -----
    For reducing at each 'near-prime' r = 22 is used as a switch between incremental reduction
    and exponential reduction, as there is a 20% chance that a number with p**22 also contains
    p**39 for n <= 499**17, and the multiplicity algorithm is most efficient when finding prime
    powers >= 17. 499**17 is not significant other than the largest prime <= 500 and 17 being
    the power at which multiplicity becomes most efficient. This value of 22 was chosen largely
    because a number was required and 22 seemed decent.
    """

    # remove as many powers of 2 as possible
    t = trailing(n)
    if t:
        n >>= t
        factors[2] = t

    r = 0
    while not n % 3:
        n //= 3
        r += 1
        if r == 22:
            r += multiplicity(3, n)
            break

    if r:
        factors[3] = r

    # similarly reduce powers of 'primes' ascending until limit
    p = 5
    while 1:

        if p > limit:
            break

        r = 0
        while not n % p:
            n //= p
            r += 1
            if r == 22:
                rr = multiplicity(p, n)
                n //= p ** rr
                r += rr
                break

        if r:
            factors[p] = r

        p += 2  # 6k + 1

        # since all smaller factors have been removed, p is a factor iff kp | n w/ k >= p
        if p * p > n:
            break

        if p > limit:
            break

        r = 0
        while not n % p:
            n //= p
            r += 1
            if r == 22:
                rr = multiplicity(p, n)
                n //= p ** rr
                r += rr
                break

        if r:
            factors[p] = r

        p += 4  # 6k - 1

        if p * p > n:
            break

    return n, p


def _check_factored(factors, n):
    """
    Checks to see if n has been completely factored, raising StopIteration if so.

    Parameter factors is directly modified and it does not need to be returned.
    """
    if isprime(n):
        if n in factors:
            factors[n] += 1
        else:
            factors[n] = 1
        n = 1

    if n == 1:
        raise StopIteration
