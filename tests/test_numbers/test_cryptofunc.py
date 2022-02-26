from cryptography318.number.crypto_functions import *
from cryptography318.number.prime import *
from cryptography318.deprecated.linear_algebra import *
from cryptography318.number.qs_bruteforce import *
from timeit import timeit
import numpy
import pytest
from sympy import factorint


@pytest.mark.skip
def test_kernel():
    from math import e

    n = randprime(10000) + 2
    L = pow(e, sqrt(log(n) * log(log(n))))
    B = int(pow(L, 1 / sqrt(2)))

    primes = primes(B)
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

        assert sqrt_mod(a, p) is not None


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
        primes = primes_gen(B := randrange(20, 50))
        result = [0]
        for p in primes:
            assert isprime(p) and result[-1] < p <= B
            result.append(p)


def test_vigenere():
    key = 'dog'
    plaintext = 'attack'
    res = vigenere_encrypt(key, plaintext)
    assert res == 'dhzavd'


def test_rho_factor():
    mix1 = lambda e: (pow(e, 2, n) + 1) % n
    mix2 = lambda e: pow(e, 2, n) + 1

    time_1, time_2 = 0, 0
    for _ in range(100):
        n = randrange(pow(10, 4), pow(10, 6))
        time_1 += timeit(lambda: pollard_rho_factor(n, mix1), number=100)
        time_2 += timeit(lambda: pollard_rho_factor(n, mix2), number=100)

    print(f"factor 1: {time_1:.2f}s")
    print(f"factor 2: {time_2:.2f}s")


def test_factor():
    time_1, time_2, time_3 = 0, 0, 0
    try:
        for _ in range(1000):
            n = randrange(pow(2, 15))
            time_1 += timeit(lambda: factor(n), number=100)
            time_2 += timeit(lambda: factor2(n), number=100)
            time_3 += timeit(lambda: factorint(n), number=100)
    except RecursionError:
        print(n)
    print(f"factor 1: {time_1:.2f}s")
    print(f"factor 2: {time_2:.2f}s")
    print(f"factor 3: {time_3:.2f}s")

    print(f"sy.factorint = 100%; factor1 = {(time_1 / time_3):.2f}%; factor2 = {(time_2 / time_3):.2f}%")


if __name__ == '__main__':
    test_factor()
