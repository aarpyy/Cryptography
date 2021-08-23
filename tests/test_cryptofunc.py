from cryptography318.crypto_functions import *
from cryptography318.tools import *
from cryptography318.prime import *
from cryptography318.linear_algebra import *
from cryptography318.quadratic_sieve import *
import timeit
import numpy
import pytest


@pytest.mark.skip
def test_kernel():
    from math import e

    n = randprime(10000) + 2
    L = pow(e, sqrt(log(n) * log(log(n))))
    B = int(pow(L, 1 / sqrt(2)))

    primes = primes_lt(B)
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


def test_lenstra():
    factors = lenstra_elliptic(55)
    print(factors)


def test_find_roots():
    for _ in range(50):
        p = randprime(pow(10, 3), pow(10, 4))
        a = randrange(2, p)
        while not quadratic_residue(a, p):
            a = randrange(2, p)

        assert find_roots(a, p) is not None


@pytest.mark.skip
def test_find_non_residue():

    # fastest method to find quadratic non-residue mod ps
    def method1():
        p = randprime(pow(10, 4), pow(10, 5))
        z = randrange(2, p)
        while not quadratic_non_residue(z, p):
            z = randrange(2, p)

    def method2():
        p = randprime(pow(10, 4), pow(10, 5))
        z = 2
        while not quadratic_non_residue(z, p):
            z += 1

    m1 = timeit.timeit(lambda: method1(), number=1000000)
    m2 = timeit.timeit(lambda: method2(), number=1000000)
    print(f"random method took {m1:.2f}s")
    print(f"increment method took {m2:.2f}s")


def test_primes_gen():
    for _ in range(50):
        primes = primes_lt_gen(B := randrange(20, 50))
        result = [0]
        for p in primes:
            assert is_prime(p) and result[-1] < p <= B
            result.append(p)


def test_vigenere():
    key = 'dog'
    plaintext = 'attack'
    res = vigenere_encrypt(key, plaintext)
    assert res == 'dhzavd'


if __name__ == '__main__':
    test_vigenere()
    caeser_shift('joxxurbqlxwbrbcnwlhrbcqnqxkpxkurwxourccunvrwmb')
