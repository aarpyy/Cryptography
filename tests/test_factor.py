from cryptography318.crypto_functions import __factor_perfect_square, pollard_rho_dlp
from cryptography318.prime import *
from cryptography318.crypto_functions import *
from cryptography318.linear_algebra import *
from cryptography318.quadratic_sieve import *
from math import prod
import time
import timeit
import pytest
import numpy


@pytest.mark.skip
def test_factor_perfect_square():
    factors = __factor_perfect_square(198103, B=20)
    n = 1
    for f in factors:
        n *= pow(f, factors[f])
    assert n == 198103


@pytest.mark.skip
def test_factor_int(power=1):
    p = randprime(pow(10, power), pow(10, power + 1)) * randprime(pow(10, power + 1), pow(10, power + 2))
    print(timeit.timeit(lambda: quadratic_sieve(p), number=100000))
    print(timeit.timeit(lambda: factor(p), number=100000))


def test_pollard_rho(it=50):
    for _ in range(it):
        g = 4
        p = randprime(pow(2, 30))
        e = randrange(p - 1)
        h = pow(g, e, p)
        r = pollard_rho_dlp(g, h, p)
        assert pow(g, r, p) == h


def test_factor_if_smooth():
    n = 5 * 7 * 7 * pow(3, 4)
    assert factor_if_smooth(n, [2, 3, 5, 7]) == [0, 4, 1, 2]
    assert factor_if_smooth(n, [2, 3, 5]) is None


def test_qs(it=10):
    for _ in range(it):
        a = randprime(pow(10, 4), pow(10, 5))
        b = randprime(pow(10, 4), pow(10, 5))
        n = a * b
        factors = quadratic_sieve(n)
        if factors is None:
            print(a, b, n)
        result = 1
        for f in factors:
            result *= pow(f, factors[f])
        assert result == n


@pytest.mark.skip
def test_max_qs():

    # this test shouldn't be run with pytest, it just purely for finding largest integer than can be factored
    # in a short amount of time by quadratic sieve
    start = time.time()
    a = randprime(pow(10, 9))
    b = randprime(pow(10, 5))
    n = a * pow(b, 2)
    print(a, b, n)
    print(quadratic_sieve(n))
    print(f"This took {time.time() - start:.2f}s")


def test_b_smooth():
    for _ in range(50):
        B = 25
        primes = primes_lt(B)
        n = 1
        for i in range(len(primes)):
            n *= pow(primes[i], randrange(2, 5))

        assert b_smooth(n, B)
        assert b_smooth(n, factors=primes)


if __name__ == '__main__':
    test_b_smooth()
