from cryptography318.prime import *
import pytest


def test_primes_lt():
    primes = primes_lt(11)
    assert primes == [2, 3, 5, 7, 11]
    assert len(primes) == prime_pi(11)

    primes = primes_lt_gen(11)
    n = 2
    for p in primes:
        assert p == n
        n = next_prime(n)


def test_is_prime():
    n = 2047
    assert not is_prime(n)

    n = randprime(pow(10, 3))
    assert confirm_prime(n)


if __name__ == '__main__':
    test_primes_lt()
    test_is_prime()
