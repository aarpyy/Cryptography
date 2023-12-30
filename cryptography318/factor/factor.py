from math import gcd, log, sqrt, log2, ceil, isqrt
from pathlib import Path
from random import randrange

from .qs import qs
from .elliptic import ecm
from cryptography318.prime import isprime, next_prime, primesieve
from cryptography318.utils.root import integer_nth_root, is_square
from cryptography318.utils.misc import as_int


def factor(n, use_rho=True, use_ecm=True, use_pm1=True, use_siqs=True, limit=None, verbose=False, *, details=None):
    """
    Attempts to factor given integer with four methods, returning None if un-factorable.
    Function first checks if number is prime then finds all small factors if any exist.
    Then (assuming no specific methods were set to False) Pollard's Rho, Lenstra's ECM,
    Pollard's P-1, and the Quadratic Sieve are used sequentially to find non-trivial factor
    of n. If any of these methods succeed, a recursive call of factor is used to factor
    the remaining n. All boolean values for use of methods are preserved in the
    recursive call.

    :param n: int number to be factored
    :param use_rho: bool determining if Pollard's Rho algorithm should be used
    :param use_ecm: bool determining if Lenstra's ECM algorithm should be used
    :param use_pm1: bool determining if Pollard's P-1 algorithm should be used
    :param use_siqs: bool determining if Quadratic Sieve algorithm should be used
    :param limit: integer limit of factors to be found using small_factors()
    :param verbose: bool determining if factorization details should be printed
    :param details: bool determining if factor details should be updated
    :return: dictionary of all primes factors and their powers, or None if not factorable
    """

    if details is None:
        details = {}
    details['methods'] = details.get('methods', [])

    n = as_int(n)

    factors = {}
    if n < 0:
        n = -n
        factors[-1] = 1

    if isprime(n):
        factors[n] = 1
        return factors

    k, _ = factor_small(factors, n, limit or 32768)

    # Update small factors in details only if we found any
    if k != n:
        details['methods'].append({
            'function': factor_small.__name__,
            'name': 'Trial division',
            'value': factors.copy()
        })

    # If we factored it completely with small factors, return factors
    if k == 1:
        return factors

    # If we factored it partially with small factors, recursively factor the rest
    if k != n:
        # If remaining factor is prime, we are done factoring
        if isprime(k):
            factors[k] = factors.get(k, 0) + 1
            return factors
        n = k

    if is_square(n):
        root = isqrt(n)
        factors[root] = 2
        details['methods'].append({
            'function': is_square.__name__,
            'name': 'Perfect square',
            'value': {root: 2}
        })
        return factors

    # Now let's check for pefect powers
    max_power = ceil(log2(n))  # Largest possible power of 2 (smallest prime) so all other powers would be smaller

    # Start loop at 3 since we can use is_square() for 2
    for i in range(3, max_power + 1):
        root = integer_nth_root(n, i)
        if root ** i == n:
            factors[root] = i
            details['methods'].append({
                'function': integer_nth_root.__name__,
                'name': 'Perfect power',
                'value': {root: i}
            })
            return factors

    factor_kwargs = {"use_rho": use_rho, "use_ecm": use_ecm, "use_pm1": use_pm1, "use_siqs": use_siqs, "limit": limit,
                     "details": details}

    if use_rho:
        value = pollard_rho_factor(n)
        if value is not None:
            details['methods'].append({
                'function': pollard_rho_factor.__name__,
                'name': "Pollard's Rho",
                'value': value
            })
            n = _factor_further(n, value, factors, **factor_kwargs)
            if n == 1:
                return factors

    if use_ecm:
        value = ecm(n, verbose=verbose)
        if value is not None:
            details['methods'].append({
                'function': ecm.__name__,
                'name': "ECM",
                'value': value
            })
            n = _factor_further(n, value, factors, **factor_kwargs)
            if n == 1:
                return factors

    if use_pm1:
        value = pollard_pm1(n)
        if value is not None:
            details['methods'].append({
                'function': pollard_pm1.__name__,
                'name': "Pollard's P-1",
                'value': value
            })
            n = _factor_further(n, value, factors, **factor_kwargs)
            if n == 1:
                return factors

    if use_siqs:
        # Use local primes.txt if we can find it, otherwise don't use a file (slow)
        path = Path(__file__).parent.parent.joinpath('data/primes.txt')
        fp = path if path.is_file() else None

        value = qs(n, fp=fp, verbose=False)
        if value is not None:
            details['methods'].append({
                'function': qs.__name__,
                'name': "Quadratic Sieve",
                'value': value
            })

            # Nothing left after quadratic sieve, so just return factors
            n = _factor_further(n, value, factors, **factor_kwargs)
            if n != 1:
                factors[n] = factors.get(n, 0) + 1
            return factors

    return factors


def pollard_pm1(n, B=None, retry=5):
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

    a = 2
    for _ in range(retry):
        m = a
        for j in primesieve.primerange(B + 1):
            exp = int(log(B, j))
            m = pow(m, pow(j, exp), n)
        q = gcd(m - 1, n)
        if 1 < q < n:
            return q

        a = randrange(2, n - 2)

    return None


def pollard_rho_factor(n, mix=None, retry=5):
    if n < 10:
        return factor_small({}, n, 10)

    if not callable(mix):
        def mix(e):
            return (pow(e, 2, n) + 1) % n

    y = 2
    for _ in range(retry):
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
