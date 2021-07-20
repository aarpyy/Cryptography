import math
import random
import sympy

from cryptography318.bailliepsw_helper import LucasPseudoPrime, D_chooser


def KnownPrime(n):
    """Helper function, confirming prime candidate is not easily known"""

    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
        if n == p:
            return True
        elif n % p == 0:
            return False
    return None


def IsPrime(n):
    """General purpose primality function that returns False if guaranteed composite"""

    # if prime candidate is too large, just performs a more accurate Miller-Rabin test
    # instead of both Baillie-PSW and Miller-Rabin;
    # returns just Baillie-PSW if less than 2^64, since Baillie-PSW is proven to be
    # deterministic for pseudo-primes up to 2^64
    if n < pow(2, 64):
        return BailliePSW_Primality(n)
    elif n > pow(2, 1000):
        return MillerRabinPrimality(n, 40)
    return BailliePSW_Primality(n) and MillerRabinPrimality(n, 10)


def MillerRabinPrimality(n, k=40):
    """MRPrimality test reduces n - 1 to a power of 2 and an odd number, then
    tests if random a is a witness of composite-ness d times, testing with
    k random a's"""

    if KnownPrime(n) is not None:
        return KnownPrime(n)

    d = n - 1
    r = 0
    while d % 2 == 0:
        r += 1
        d >>= 1

    for _ in range(k):
        if not MillerTest(d, n):
            return False

    return True


def MillerTest(d, n):
    """Helper function for MRPrimality which uses previously found d to
    check if random a is a witness to n's composite-ness"""

    a = random.randrange(2, n - 1)
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


def MillerRabin_base_a(a, n):
    """Miller Rabin test with specific base of a"""

    if math.gcd(a, n) > 1:
        return False
    q = n - 1
    k = 0
    while q % 2 == 0:
        q >>= 1
        k += 1

    a = pow(a, q, n)
    if a == 1 or a == n - 1:
        return True
    for _ in range(k):
        if a == -1 or a == n - 1:
            return True
        elif a == 1:
            return False
        a = pow(a, 2, n)

    return False


def RandomPrime(*args):
    """Uses combination of Miller-Rabin and Baillie-PSW primality tests to generate random prime"""

    base_2 = False

    # determines if user entered a lower and upper limit or just an upper
    if len(args) not in [1, 2]:
        raise TypeError("Usage: RandomPrime(limit=int) or RandomPrime(base=int, limit=int)")
    base, limit = (args[0], args[1]) if len(args) == 2 else (3, args[0])

    if base == 2:
        base_2 = True

    base = base | 1

    # if base_2, uses 2 as a base and increments by 1 (default) for generating random int
    # if base =/= 2, generates random int starting at lower limit, incrementing by 2
    while True:
        prime = random.randrange(2, limit) if base_2 else random.randrange(base, limit, 2)
        if IsPrime(prime):
            return prime


def ConfirmPrime(n):
    """Uses infinitely deterministic primality test, checking if candidate has factors
    of any primes <= square root of candidate"""

    if KnownPrime(n) is not None:
        return KnownPrime(n)

    for num in range(3, (math.isqrt(n) + 1) | 1, 2):
        witness = list(map(lambda x: MillerRabin_base_a(x, n), [2, 3, 5, 7, 11, 13]))
        if False not in witness and n % num == 0:
            return False
    return True


def NextPrime(n):
    """Returns first prime after number given"""

    # ensures n is odd to start so that can increment by 2
    n = n | 1
    while True:
        if IsPrime(n):
            return n
        n += 2


def PrevPrime(n):
    """Returns first prime before number given"""

    # ensures n is odd to start so that can decrement by 2
    n = n | 1
    while True:
        if IsPrime(n):
            return n
        n -= 2


def BailliePSW_Primality(candidate):
    """Perform the Baillie-PSW probabilistic primality test on candidate"""

    # Check divisibility by a short list of primes less than 50
    if KnownPrime(candidate) is not None:
        return KnownPrime(candidate)

    # Now perform the Miller-Rabin primality test base 2
    if not MillerRabin_base_a(2, candidate):
        return False

    # Checks if number has integer square root, if it does not, math.isqrt
    # will not be exact square root, if it has integer square root then
    # math.isqrt will square perfectly to candidate
    if math.isqrt(candidate) ** 2 == candidate:
        return False

    # Finally perform the Lucas primality test
    D = D_chooser(candidate)
    if not LucasPseudoPrime(candidate, D, 1, (1 - D) / 4):
        return False

    return True


def PollardP1(n, limit=pow(10, 5)):
    """Pollard's p - 1 algorithm for factoring large composites.
    Returns a factor if factor-able, False if otherwise."""

    if IsPrime(n):
        raise ValueError("Make sure to enter a composite number")

    for a in [2, 3, 5]:
        m = a
        for j in range(2, limit):
            m = pow(m, j, n)
            k = math.gcd(m - 1, n)
            if 1 < k < n:
                return k

    return False
