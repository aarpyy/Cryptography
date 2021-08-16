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
    p = RandomPrime(pow(10, power), pow(10, power + 1)) * RandomPrime(pow(10, power + 1), pow(10, power + 2))
    print(timeit.timeit(lambda: quadratic_sieve(p), number=100000))
    print(timeit.timeit(lambda: FactorInt(p), number=100000))


def test_pollard_rho(it=50):
    for _ in range(it):
        g = 4
        p = RandomPrime(pow(2, 30))
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
        a = RandomPrime(pow(10, 4), pow(10, 5))
        b = RandomPrime(pow(10, 4), pow(10, 5))
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
    a = RandomPrime(pow(10, 9))
    b = RandomPrime(pow(10, 5))
    n = a * pow(b, 2)
    print(a, b, n)
    print(quadratic_sieve(n))
    print(f"This took {time.time() - start:.2f}s")


@pytest.mark.skip
def test_kernel():
    from math import e

    n = RandomPrime(10000) + 2
    L = pow(e, sqrt(log(n) * log(log(n))))
    B = int(pow(L, 1 / sqrt(2)))

    primes = PrimesLT(B)
    bases, squares, exp = find_perfect_squares(n, primes)
    matrix = Matrix(exp, mod=2).astype(numpy.int64)
    print("done with first part")

    def kern1(m):
        kernel(m)

    def kern2(m):
        kern1(m)

    original = timeit.timeit(lambda: kern1(matrix), number=100000)
    print(f"original took {original:.2f}s")
    new = timeit.timeit(lambda: kern2(matrix), number=100000)
    print(f"new took {new:.2f}s")


def test_find_perfect_squares():
    n = 100
    primes = [2, 3, 5, 7]
    base, sq, exp = find_perfect_squares(n, primes)
    for i, b in enumerate(base):
        assert exp[i] == factor_if_smooth(pow(b, 2) - n, primes)
        assert exp_value(exp[i], primes=primes) == sq[i]


@pytest.mark.skip
def test_all_tests():
    test_factor_perfect_square()
    test_factor_int()


def test_lenstra():
    factors = lenstra_elliptic(55)
    print(factors)


if __name__ == '__main__':
    test_qs()
