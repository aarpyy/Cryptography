import random, math
from crypto_functions import *


# from internet
def IsPrime(num):
    # 2 is the only even prime, checks for 2 first
    if num == 2:
        return True
    # num & 1 returns intersection of binary 1 and binary num
    # if num is odd, it will always intersect with 1 and return True
    # not num & 1 filters all evens to return False, otherwise check below
    if not num & 1:
        return False
    # checks if fermat's little thereom works, will sometimes produce psuedo-primes
    return pow(2, num - 1, num) == 1


def MillerRabinPrimality(n, k=40):
    if n == 2 or n == 3 or n == 5:
        return True
    elif n % 2 == 0:
        return False
    else:
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
    a = random.randrange(2, n - 1)

    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return True

    while d != n - 1:
        x = pow(x, 2, n)
        d <<= 1

        if x == 1:
            return False
        elif x == n - 1:
            return True
    return False


def MR(a, n):
    if GCD(a, n) > 1:
        return "Composite"
    q = n-1
    k = 0
    while q % 2 == 0:
        q //= 2
        k += 1

    a = pow(a, q, n)
    if a == 1:
        return "Test Fails"

    for _ in range(k):
        if a == -1 or a == n - 1:
            return "Test Fails"
        a = pow(a, 2, n)

    return "Composite"


# certainty value represents probability; if k = certainty value,
# probability that number generated is prime = 4 ^ -k
def RandomPrime(base, limit=None, certainty=40):
    base_2 = False

    # determines if user entered a lower and upper limit or just an upper
    if limit is not None:
        if base == 2:
            base_2 = True
        base = base | 1
    else:
        limit = base
        base = 3

    # if base_2, uses 2 as a base and increments by 1 (default) for generating random int
    if base_2:
        while True:
            prime = random.randrange(2, limit)
            if MillerRabinPrimality(prime, certainty):
                return prime

    # if base =/= 2, generates random int starting at lower limit, incrementing by 2
    while True:
        prime = random.randrange(base, limit, 2)
        if MillerRabinPrimality(prime, certainty):
            return prime


def ConfirmPrime(n):
    if n == 2:
        return True
    elif not n & 1:
        return False
    elif n % 3 == 0 or n % 5 == 0:
        return False
    else:
        sq = math.floor(math.sqrt(n))
        primes = []
        for num in range(3, sq):
            if IsPrime(num):
                primes.append(num)

        for num in primes:
            if n % num == 0:
                return False
        return True


def NextPrime(n):
    n = n | 1

    while True:
        if IsPrime(n):
            if MillerRabinPrimality(n):
                return n
        n += 2


# https://rosettacode.org/wiki/Jacobi_symbol#Python - modified
def Jacobi(a, n):
    assert n > 0, n & 1
    a %= n
    result = 1
    while a != 0:
        while a % 2 == 0:
            a >>= 1
            n_mod_8 = n % 8
            if n_mod_8 in (3, 5):
                result *= -1
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result *= -1
        a %= n
    if n == 1:
        return result
    else:
        return 0