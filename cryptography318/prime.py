from math import gcd, isqrt
from random import randrange
import ctypes
import pathlib

from cryptography318.bailliepsw_helper import LucasPseudoPrime, D_chooser


def KnownPrime(n, first_n=50):
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

    if not 2 < first_n < 167:
        first_n = known_primes
    else:
        first_n = known_primes[:first_n]

    for p in first_n:
        if n == p:
            return True
        elif n % p == 0:
            return False
    return None


def IsPrime(n):
    """
    IsPrime function returns False iff the prime-candidate is composite, and True
    if the prime-candidate is probably prime.

    Uses deterministic variants of the Miller-Rabin Primality test, which, through
    the use of specific bases and ranges, can deterministically return True iff
    candidate is prime for n < 3317044064679887385961981. For all larger n,
    there is no  known set of bases that makes the MR test deterministic. Thus a
    SPRP-test consisting of a Strong Lucas Pseudo-prime test and a Miller-Rabin
    test with 20 random bases a, s.t. 1 < a < n is used to determine if candidate is
    probably prime.
    """

    if KnownPrime(n, first_n=20) is not None:
        return KnownPrime(n)

    if n < 2047:
        return MillerRabin_bases([2], n)
    if n < 1373653:
        return MillerRabin_bases([2, 3], n)
    if n < 9080191:
        return MillerRabin_bases([31, 73], n)
    if n < 1050535501:
        return MillerRabin_bases([336781006125, 9639812373923155], n)
    if n < 3215031751:
        return MillerRabin_bases([2, 3, 5, 7], n)
    if n < 4759123141:
        return MillerRabin_bases([2, 7, 61], n)
    if n < 1122004669633:
        return MillerRabin_bases([2, 13, 23, 1662803], n)
    if n < 55245642489451:
        return MillerRabin_bases([2, 141889084524735, 1199124725622454117, 11096072698276303650], n)
    if n < 7999252175582851:
        return MillerRabin_bases([2, 4130806001517, 149795463772692060, 186635894390467037, 3967304179347715805], n)
    if n < 18446744073709551616:
        return MillerRabin_bases([2, 325, 9375, 28178, 450775, 9780504, 1795265022], n)
    if n < 318665857834031151167461:
        return MillerRabin_bases([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37], n)
    if n < 3317044064679887385961981:
        return MillerRabin_bases([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41], n)

    return MillerRabinPrimality(n, k=40) and BailliePSW_Primality(n, mr=False)


def MillerRabinPrimality(n, k=40):
    """MRPrimality test reduces n - 1 to a power of 2 and an odd number, then
    tests if random a is a witness of n's composite-ness, testing with
    k random a's"""

    if res := KnownPrime(n) is not None:
        return res

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


def MillerRabin_bases(lst_bases, n):
    """Helper function that allows for a list of witnesses to be tested
    using MillerRabin_base_a function"""

    for a in lst_bases:
        if not MillerRabin_base_a(a, n):
            return False
    return True


def MillerRabin_base_a(a, n):
    """Miller Rabin test with specific base of a"""

    if a >= n:
        a %= n

    if a == 0:
        return True

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
        raise ValueError("Usage: RandomPrime(limit->int) or RandomPrime(base->int, limit->int)")
    base, limit = (args[0], args[1]) if len(args) == 2 else (3, args[0])

    if base == 2:
        base_2 = True

    base = base | 1

    # if base_2, uses 2 as a base and increments by 1 (default) for generating random int
    # if base =/= 2, generates random int starting at lower limit, incrementing by 2
    while True:
        prime = randrange(2, limit) if base_2 else randrange(base, limit, 2)
        if IsPrime(prime):
            return prime


def AllFactors(n):
    """Uses infinitely deterministic primality test, checking if candidate has factors
    of any primes <= square root of candidate"""

    if res := KnownPrime(n) is not None:
        return res

    for num in range(3, (isqrt(n) + 1) | 1, 2):
        witness = list(map(lambda x: MillerRabin_base_a(x, n), [2, 3, 5, 7, 11, 13]))
        if False not in witness and n % num == 0:
            return False
    return True


def ConfirmPrime(n):
    """Uses infinitely deterministic AKS (Agrawal-Kayal-Saxena) primality test which
    returns True if-and-only-if n is prime"""

    if res := KnownPrime(n) is not None:
        return res

    # generates the n-th row of Pascal's triangle, if any of the coefficients != 0 mod n, n is not prime
    for k in range(1, (n+1)//2):
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
    n = (n + 1) | 1
    while True:
        if IsPrime(n):
            return n
        n += 2


def PrevPrime(n):
    """Returns first prime before number given"""

    # ensures n is odd to start so that can decrement by 2
    n = (n - 2) | 1
    while True:
        if IsPrime(n):
            return n
        n -= 2


def BailliePSW_Primality(candidate, mr=True):
    """Perform the Baillie-PSW probabilistic primality test on candidate"""

    # Check divisibility by a short list of primes less than 50
    if KnownPrime(candidate) is not None:
        return KnownPrime(candidate)

    # Now perform the Miller-Rabin primality test base 2
    if mr and not MillerRabin_base_a(2, candidate):
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
