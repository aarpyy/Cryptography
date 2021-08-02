from cryptography318.crypto_functions import _factorPerfectSquare
from cryptography318.prime import *
from cryptography318.crypto_functions import *
from math import prod


def test_factor_perfect_square():
    factors = _factorPerfectSquare(198103, B=20)
    n = 1
    for f in factors:
        n *= pow(f, factors[f])
    assert n == 198103


def test_factor_int(power=6):
    p = RandomPrime(pow(10, power), pow(10, power + 1)) * RandomPrime(pow(10, power + 1), pow(10, power + 2))
    factors = FactorInt(p)
    assert p == prod(map(lambda b, e: pow(b, e), factors.keys(), factors.values()))


def test_all_tests():
    test_factor_perfect_square()
    test_factor_int()


if __name__ == '__main__':
    test_all_tests()
