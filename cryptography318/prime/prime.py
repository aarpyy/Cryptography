from math import prod
from random import randrange, randint

from cryptography318.prime.bailliepsw_helper import D_chooser, LucasPseudoPrime
from cryptography318.utils.root import is_square
from cryptography318.prime.primesieve import primesieve


def miller_rabin(n, k=40, *, details=False):
    """
    MRPrimality test reduces n - 1 to a power of 2 and an odd number, then
    tests if random `a` is a witness of n's composite-ness, testing with
    k random a's
    """

    d = n - 1
    r = 0
    while not d & 1:
        r += 1
        d >>= 1

    for _ in range(k):
        if not _mr_test(d, n):
            details['methods'].append({
                'function': _mr_test.__name__,
                'name': 'Miller-Rabin',
                'description': f"{d} is a witness to {n}'s composite-ness",
                'value': False
            })
            return False

    details['methods'].append({
        'function': _mr_test.__name__,
        'name': 'Miller-Rabin',
        'description': f"Using {k} random bases, {n} is probably prime",
        'value': True
    })
    return True


def _mr_test(d, n):
    """Helper function for miller_rabin which uses previously found d to
    check if random `a` is a witness to n's composite-ness"""

    a = randrange(2, n - 1)
    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return True

    # doubles d every time until d returns to original n-1 value
    while d != n - 1:
        x = pow(x, 2, n)
        d <<= 1

        if x == 1:
            return False
        elif x == n - 1:
            return True
    return False


def _miller_rabin_base_a(a, n):
    """Miller Rabin test with specific base of a"""

    if a >= n:
        a %= n

    if not a:
        return True

    q = n - 1
    k = 0
    while not q & 1:
        q >>= 1
        k += 1

    a = pow(a, q, n)
    if a == 1 or a == n - 1:
        return True
    for _ in range(k):
        # If we found any a^2 = -1 mod n then we know `a` is not a witness to n's compositeness
        if a == -1 or a == n - 1:
            return True

        # If we found an a^2 = 1 mod n where a != +/- 1 then a is definitely composite
        elif a == 1:
            return False
        a = pow(a, 2, n)

    return False


def miller_rabin_bases(bases, n, *, details=None):
    """Helper function that allows for a list of witnesses to be tested
    using MillerRabin_base_a function"""

    if details is None:
        details = {}
    details['methods'] = details.get('methods', [])

    for a in bases:
        if not _miller_rabin_base_a(a, n):
            details['methods'].append({
                'function': _miller_rabin_base_a.__name__,
                'name': 'Miller-Rabin',
                'description': f"{a} is a witness to {n}'s composite-ness",
                'value': False
            })
            return False

    details['methods'].append({
        'function': _miller_rabin_base_a.__name__,
        'name': 'Miller-Rabin',
        'description': f"Using {', '.join(str(b) for b in bases)} as bases, {n} is probably prime",
        'value': True
    })
    return True


def baillie_psw(n, mr=True, details=None):
    """
    Perform the Baillie-PSW probabilistic primality test on candidate.

    :param n: prime candidate
    :param mr: if Miller-Rabin test base 2 should be used
    :param details:
    :return:
    """
    if details is None:
        details = {}
    details['methods'] = details.get('methods', [])

    # Check divisibility by a short list of primes less than 50
    if (res := known_prime(n)) is not None:
        details['methods'].append({
            'function': known_prime.__name__,
            'name': 'Known prime',
            'value': res
        })
        return res

    # Now perform the Miller-Rabin primality test base 2
    if mr and not _miller_rabin_base_a(2, n):
        details['methods'].append({
            'function': _miller_rabin_base_a.__name__,
            'name': 'Miller-Rabin',
            'description': f"2 is a witness to {n}'s composite-ness",
            'value': False
        })
        return False

    # Checks if number has square root
    if is_square(n):
        details['methods'].append({
            'function': is_square.__name__,
            'name': 'Is square',
            'value': False
        })
        return False

    # Finally perform the Lucas primality test
    D = D_chooser(n)
    value = LucasPseudoPrime(n, D, 1, (1 - D) // 4)
    details['methods'].append({
        'function': LucasPseudoPrime.__name__,
        'name': 'Lucas pseudo-prime',
        'value': value
    })
    return value


def known_prime(n):
    """Helper function, confirming prime candidate is not easily known"""

    known_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
                    103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,
                    211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317,
                    331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443,
                    449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577,
                    587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701,
                    709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839,
                    853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983,
                    991, 997]

    for p in known_primes:
        if n == p:
            return True
        elif n % p == 0:
            return False
    return None


def isprime(n, *, details=None):
    """
    IsPrime function returns False iff the prime-candidate is composite, and True
    if the prime-candidate is probably prime.

    Uses deterministic variants of the Miller-Rabin Primality test, which, through
    the use of specific bases and ranges, can deterministically return True iff
    candidate is prime for n < 3317044064679887385961981. For all larger n,
    there is no  known set of bases that makes the MR test deterministic. Thus, a
    SPRP-test consisting of a Strong Lucas Pseudo-prime test and a Miller-Rabin
    test with 20 random bases `a`, s.t. 1 < a < n is used to determine if candidate is
    probably prime.
    """

    if details is None:
        details = {}
    details['methods'] = details.get('methods', [])

    if n < 2:
        details['error'] = str(ValueError("n must be greater than 1"))
        return False

    elif n < 10:
        is_prime = bool([0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0][n])
        details['methods'].append({
            'function': None,
            'name': 'Known prime' if is_prime else 'Known composite',
            'value': is_prime
        })
        return is_prime

    # check for odds
    elif not n & 1:
        details['methods'].append({
            'function': None,
            'name': 'Even',
            'value': False
        })
        return False

    # check for all other instances n != 6k +/- 1
    elif not n % 3:
        details['methods'].append({
            'function': None,
            'name': 'Divisible by 3',
            'description': 'All prime integers are of the form 6k +/- 1 for some positive integer k',
            'value': False
        })
        return False

    # This step is pretty useless unless primesieve is being used for something else or is
    # being purposefully generated, since it is constructed only with first 6 primes
    if n < primesieve.tail and n in primesieve:
        details['methods'].append({
            'function': None,
            'name': 'Cached prime',
            'value': True
        })
        return True
    elif n < 2047:
        return miller_rabin_bases([2], n, details=details)
    elif n < 1373653:
        return miller_rabin_bases([2, 3], n, details=details)
    elif n < 9080191:
        return miller_rabin_bases([31, 73], n, details=details)
    elif n < 1050535501:
        return miller_rabin_bases([336781006125, 9639812373923155], n, details=details)
    elif n < 3215031751:
        return miller_rabin_bases([2, 3, 5, 7], n, details=details)
    elif n < 4759123141:
        return miller_rabin_bases([2, 7, 61], n, details=details)
    elif n < 1122004669633:
        return miller_rabin_bases([2, 13, 23, 1662803], n, details=details)
    elif n < 55245642489451:
        return miller_rabin_bases([2, 141889084524735, 1199124725622454117, 11096072698276303650], n, details=details)
    elif n < 7999252175582851:
        return miller_rabin_bases([2, 4130806001517, 149795463772692060, 186635894390467037, 3967304179347715805], n,
                                  details=details)
    elif n < 18446744073709551616:
        return miller_rabin_bases([2, 325, 9375, 28178, 450775, 9780504, 1795265022], n, details=details)
    elif n < 318665857834031151167461:
        return miller_rabin_bases([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37], n, details=details)
    elif n < 3317044064679887385961981:
        return miller_rabin_bases([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41], n, details=details)
    else:
        res = miller_rabin(n, k=40, details=details)
        if not res:
            return False

        value = baillie_psw(n, mr=False, details=details)
        # If miller rabin didn't say it was composite, then we should note that we used Baillie-PSW
        details['methods'].append({
            'function': baillie_psw.__name__,
            'name': 'Baillie-PSW',
            'description': 'Baillie-PSW is a combination of Miller-Rabin and Lucas Pseudo-Prime tests',
            'value': value
        })
        return value


def randprime(a: int, b: int = None):
    """Uses combination of Miller-Rabin and Baillie-PSW primality tests to generate random prime

    Note
    ----
    If no lower bound is specified, 2 will never be generated, since bounds of sampling are set [3, b)

    :param a: integer starting point of range for random prime
    :param b: integer stopping point of range for random prime (exclusive)
    """

    # determines if user entered a lower and upper limit or just an upper
    if b is None:
        a, b = 2, a

    if a >= b:
        return
    a, b = map(int, (a, b))
    n = randint(a - 1, b)
    p = next_prime(n)
    if p >= b:
        p = prev_prime(b)
    if p < a:
        raise ValueError("no primes exist in the specified range")
    return p


def confirm_prime(n):
    """Uses deterministic AKS (Agrawal-Kayal-Saxena) primality test which
    returns True if-and-only-if n is prime"""

    if n < 2:
        return False
    elif n == 2:
        return True

    if n in primesieve:
        return True

    # generates the n-th row of Pascal's triangle, if any of the coefficients != 0 mod n, n is not prime
    for k in range(1, (n + 1) // 2):
        res = 1
        if k > (n - k):
            k = n - k
        for i in range(k):
            res *= n - i
            res //= i + 1
        if res % n != 0:
            return False
    return True


def next_prime(n):
    """Returns first probable prime after number given"""

    if n < 2:
        return 2
    elif n < 11:
        return [2, 2, 3, 5, 5, 7, 7, 11, 11, 11, 11][n]
    elif n < primesieve.tail:
        a, b = primesieve.search(n)
        return primesieve[b]

    # Ensures that n starts at the nearest 6k + 1

    # `a` is the closest 6k + 1 to n
    a = n - (n % 6) + 1
    if a <= n:

        # If a <= n, try 6k + 5, only return if that's greater than n
        a += 4
        if a > n and isprime(a):
            return a

        # Otherwise, get `a` to 6k + 1, if this is prime, it's guaranteed > n so return
        a += 2
        if isprime(a):
            return a

        # Otherwise, since 6k + 1 above n didn't work, set n to 6k + 5 for below loop
        n = a + 4
    elif isprime(a):

        # If a > n and is prime, just return (this case only runs when n % 6 == 0 and n + 1 is prime)
        return a
    else:

        # Otherwise, start off n at a + 4 which is (6k + 1) + 4
        n = a + 4

    assert n % 6 == 5

    # Iterate up through each 6k +/- 1
    while True:
        if isprime(n):
            return n
        n += 2
        if isprime(n):
            return n
        n += 4


def prev_prime(n):
    """Returns first prime before number given"""

    if n < 3:
        raise ValueError(f"No primes exist < {n}")

    if n < 11:
        return [0, 0, 0, 0, 3, 3, 5, 5, 7, 7, 7][n]

    if n <= primesieve.tail:
        i = primesieve.search(n)
        if isinstance(i, tuple):
            i = i[0]
        else:
            i -= 1
        return primesieve[i]

    # ensures that n starts at the nearest 6k - 1 below
    r = n % 6
    if not r:
        n -= 1
    elif r == 1:
        n -= 2
    elif r <= 4:
        n -= r
        if isprime(n + 1):
            return n + 1
        n -= 1
    elif r == 5:
        n -= 4
        if isprime(n):
            return n
        n -= 2

    while True:
        if isprime(n):
            return n
        n -= 4
        if isprime(n):
            return n
        n -= 2

        # If we are sure that n passed was larger than 2, and we somehow end up below 2, just return 2
        if n < 2:
            return 2


def prime_range(a, b=None):
    """Constructs list of a <= primes < b"""
    if b is None:
        b = a
        a = 1

    if b < 2:
        return []

    global primesieve
    if b < primesieve.tail:
        primesieve.extend(b)

    if a < 2:
        a = 2
    indices = primesieve.search(a, b)
    start, stop = indices[0], indices[1]
    if isinstance(start, tuple):
        start = start[1]
    if isinstance(stop, tuple):
        stop = stop[1]
    return primesieve[start:stop]


def sqrt_mod(a, p):
    """Finds a solution for x to equation x^2 = a (mod p). If a solution is returned, a second
    solution s2 will also exist where s2 = -x (mod p)."""

    if not a:
        return 0
    elif not quadratic_residue(a, p):
        return None

    mod8 = p % 8
    if mod8 == 1:
        q = p - 1
        s = 0
        while not q & 1:
            q >>= 1
            s += 1

        z = randrange(2, p)
        while not quadratic_non_residue(z, p):
            z = randrange(2, p)

        m = s
        c = pow(z, q, p)
        t = pow(a, q, p)
        r = pow(a, (q + 1) // 2, p)

        while True:
            if t == 0:
                return 0
            if t == 1:
                return r

            i = 0
            x = t
            while x != 1:
                x = pow(x, 2, p)
                i += 1

            b = pow(c, pow(2, m - i - 1), p)
            c = pow(b, 2, p)
            m = i

            t = (t * c) % p
            r = (r * b) % p
    elif mod8 == 5:
        a2 = a + a
        v = pow(a2, (p - 5) // 8, p)
        i = (a2 * v * v) % p
        return (a * v * (i - 1)) % p
    else:
        return pow(a, (p + 1) // 4, p)


def lift_sqrt(root, n, modulus, q=None):
    """
    Given integer root that is the modular square root of ``n
    % modulus`` compute and return the square root of ``n % modulus * q``.
    That is, if q is not given, this function computes the square root
    of ``n % modulus ** 2``.

    Note
    ----
    If q is not given, modulus must be prime. If q is given it must be
    prime and modulus must be a prime power of q.

    References
    ----------
    This algorithm is Hensel's Lemma.

    :param root: square root of n mod modulus
    :param n: integer in FF(modulus)
    :param modulus: modulus of field
    :param q: prime
    :return: square root of n mod modulus * q
    """
    if q is None:
        q = modulus

    s = ((n - root * root) // modulus) * pow(root + root, -1, q)
    return (root + s * modulus) % (modulus * q)


def quadratic_residue(a, p):
    """Returns True if n is a quadratic residue mod p, False otherwise. Uses Euler's criterion to assess values.
    Assumes p is odd prime."""

    return pow(a, (p - 1) // 2, p) == 1


def quadratic_non_residue(a, p):
    """Returns True if n is a quadratic non-residue mod p, False otherwise. Uses Euler's criterion to assess values.
    Assumes p is odd prime."""

    return pow(a, (p - 1) // 2, p) == p - 1


def chinese_remainder(values, moduli):
    # Initializes lists of moduli, mod = product of all moduli
    mod = prod(moduli)

    # Maps list of moduli and their inverses to x and y respectively
    x, y = [], []
    for m in moduli:
        mi = mod // m
        x.append(mi)
        y.append(pow(mi, -1, m))

    # Accumulates product of number and moduli and their inverses
    acc = 0
    for i in range(len(values)):
        acc = (acc + values[i] * x[i] * y[i]) % mod

    return acc
