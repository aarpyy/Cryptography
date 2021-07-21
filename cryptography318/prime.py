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
    if n < 2047:
        return MillerRabin_bases([2], n)
    if n < 1373653:
        return MillerRabin_bases([2, 3], n)
    if n < 9080191:
        return MillerRabin_bases([31, 73], n)
    if n < 25326001:
        return MillerRabin_bases([2, 3, 5], n)
    if n < 3215031751:
        return MillerRabin_bases([2, 3, 5, 7], n)
    if n < 4759123141:
        return MillerRabin_bases([2, 7, 61], n)
    if n < 1122004669633:
        return MillerRabin_bases([2, 13, 23, 1662803], n)
    if n < 2152302898747:
        return MillerRabin_bases([2, 3, 5, 7, 11], n)
    if n < 3474749660383:
        return MillerRabin_bases([2, 3, 5, 7, 11, 13], n)
    if n < 341550071728321:
        return MillerRabin_bases([2, 3, 5, 7, 11, 13, 17], n)
    if n < 3825123056546413051:
        return MillerRabin_bases([2, 3, 5, 7, 11, 13, 17, 19, 23], n)
    if n < 18446744073709551616:
        return MillerRabin_bases([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37], n)
    if n < 318665857834031151167461:
        return MillerRabin_bases([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37], n)
    if n < 3317044064679887385961981:
        return MillerRabin_bases([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41], n)

    return MillerRabinPrimality(n, 20) and BailliePSW_Primality(n)


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


def MillerRabin_bases(lst, n):
    """Helper function that allows for a list of witnesses to be tested
    using MillerRabin_base_a function"""

    for a in lst:
        if not MillerRabin_base_a(a, n):
            return False
    return True


def MillerRabin_base_a(a, n):
    """Miller Rabin test with specific base of a"""

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


def AllFactors(n):
    """Uses infinitely deterministic primality test, checking if candidate has factors
    of any primes <= square root of candidate"""

    if KnownPrime(n) is not None:
        return KnownPrime(n)

    for num in range(3, (math.isqrt(n) + 1) | 1, 2):
        witness = list(map(lambda x: MillerRabin_base_a(x, n), [2, 3, 5, 7, 11, 13]))
        if False not in witness and n % num == 0:
            return False
    return True


def ConfirmPrime(n):
    """Uses infinitely deterministic AKS (Agrawal-Kayal-Saxena) primality test which
    returns True if-and-only-if n is prime"""
    if KnownPrime(n) is not None:
        return KnownPrime(n)

    # generates the n-th row of Pascal's triangle, if any of the coefficients != 0 mod n, n is not prime
    for k in range(1, n):
        res = 1
        if k > (n - k):
            k = n - k
        for i in range(0, k):
            res = res * (n - i)
            res = res // (i + 1)
        if res % n != 0:
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

    # Checks if number has square root using sympy function
    from sympy.ntheory.primetest import is_square
    if is_square(candidate):
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
