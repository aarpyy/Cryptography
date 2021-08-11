from cryptography318.prime import *
import pytest


def test_primes_lt():
    primes = PrimesLT(11)
    assert primes == [2, 3, 5, 7, 11]
    assert len(primes) == PrimePi(11)

    primes = PrimesLT_gen(11)
    n = 2
    for p in primes:
        assert p == n
        n = NextPrime(n)


def test_is_prime():
    n = 2047
    assert not IsPrime(n)

    n = RandomPrime(pow(10, 3))
    assert ConfirmPrime(n)


if __name__ == '__main__':
    test_primes_lt()
    test_is_prime()
