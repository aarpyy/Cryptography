import os
from math import gcd, log, sqrt
from random import randrange

from cryptography318.factor import lenstra_ecm
from cryptography318.factor import siqs
from cryptography318.prime.prime import isprime, next_prime, primesieve


def factor(n, rho=True, ecm=True, p1=True, qs=True, limit=None, *, details=None):
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
    :param details: bool determining if factor details should be updated
    :return: dictionary of all primes factors and their powers, or None if not factorable
    """

    if details is None:
        details = {}
    details['methods'] = details.get('methods', [])

    if not isinstance(n, int):
        details['error'] = str(TypeError("n must be an integer"))
        return {}

    if n <= 1:
        details['error'] = str(ValueError("n must be greater than 1"))
        return {}

    if isprime(n):
        details['methods'].append({
            'function': isprime.__name__,
            'name': 'Is Prime',
            'value': n
        })
        return {n: 1}

    if limit is None:
        limit = 32768

    factors = {}
    k, _ = factor_small(factors, n, limit)

    # Update small factors in details - should only already exist if called recursively
    details['methods'].append({
        'function': factor_small.__name__,
        'name': 'Trial Division',
        'value': factors.copy()
    })

    # If we factored it completely with small factors, return factors
    if k == 1:
        return factors

    # If we factored it partially with small factors, recursively factor the rest
    elif k != n:
        # If remaining factor is prime, we are done factoring
        if isprime(k):
            factors[k] = factors.get(k, 0) + 1

            # Update details from previous small factors to include this prime
            details['methods'].append({
                'function': isprime.__name__,
                'name': 'Is Prime',
                'value': k
            })
            return factors
        n = k

    factor_kwargs = {"rho": rho, "ecm": ecm, "p1": p1, "qs": qs, "limit": limit, "details": details}

    if rho:
        value = pollard_rho_factor(n)
        details['methods'].append({
            'function': pollard_rho_factor.__name__,
            'name': "Pollard's Rho",
            'value': value
        })
        n = _factor_further(n, value, factors, **factor_kwargs)
        if n == 1:
            return factors

    if ecm:
        value = lenstra_ecm(n)
        details['methods'].append({
            'function': lenstra_ecm.__name__,
            'name': "ECM",
            'value': value
        })
        n = _factor_further(n, value, factors, **factor_kwargs)
        if n == 1:
            return factors

    if p1:
        value = pollard_p1(n)
        details['methods'].append({
            'function': pollard_p1.__name__,
            'name': "Pollard's P-1",
            'value': value
        })
        n = _factor_further(n, value, factors, **factor_kwargs)
        if n == 1:
            return factors

    if qs:
        # Use local primes.txt if we can find it, otherwise don't use a file (slow)
        fp = "primes.txt" if os.path.exists("../prime/primes.txt") else None

        value = siqs(n, fp=fp, loud=False)
        details['methods'].append({
            'function': siqs.__name__,
            'name': "Quadratic Sieve",
            'value': value
        })
        # Nothing left after quadratic sieve, so just return factors
        n = _factor_further(n, value, factors, **factor_kwargs)
        if n != 1:
            factors[n] = factors.get(n, 0) + 1
        return factors

    return None


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
        def mix(e):
            return (pow(e, 2, n) + 1) % n

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

        def mix(e):
            return (pow(e, 2, n) + a) % n

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
    if f <= 0:
        return n

    n //= f
    if isprime(f):
        factors[f] = factors.get(f, 0) + 1
    else:
        factors_f = factor(f, **kwargs)
        for prime in factors_f:
            factors[prime] = factors.get(prime, 0) + factors_f[prime]

    if n == 1:
        return n

    more_facs = factor(n, **kwargs)
    if more_facs is None:
        return n

    for prime in more_facs:
        factors[prime] = factors.get(prime, 0) + more_facs[prime]
    return 1
