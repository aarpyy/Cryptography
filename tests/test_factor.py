from cryptography318.crypto_functions import pollard_rho_dlp
from cryptography318.prime import *
from cryptography318.crypto_functions import *
from cryptography318.crypto_functions import _factor_with_known
from cryptography318.linear_algebra import *
from cryptography318.quadratic_sieve import *
from math import prod
import time
import timeit
import pytest
import numpy


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
    # n = randprime(pow(10, 4), pow(10, 6)) + 2
    n = 888135

    # print(f"n: {n}")
    # factors = quadratic_sieve(n, force=20)
    factors = {5: 1, 177627: 3, 16: 2}

    final = reduce(
        lambda i, c: join_dict(i, reduce(
            lambda a, b: join_dict(a, {b: k[b] * factors[c]}), k := factor(c), {}
        )) if not isprime(c) else join_dict(i, {c: factors[c]}), factors, {}
    )
    print(_factor_with_known(factors))
    # temp = reduce(lambda i, c: join_dict(i, {c: k[c] * factors[c]}), k := factor(c), {})
    # print("temp")
    # print(temp)

    more_factors = {}
    for f in factors:
        if not isprime(f):  # if further factoring to do
            more_factors[f] = factor(f)

            for e in more_factors[f]:
                more_factors[f][e] *= factors[f]  # if non-prime factor had power greater than 1, * powers of its prime factors by exp

    final = {}
    for f in factors:
        if f not in more_factors:
            final[f] = factors[f]

    final = reduce(lambda i, c: join_dict(i, more_factors[c]), more_factors, final)
    # print(final)
