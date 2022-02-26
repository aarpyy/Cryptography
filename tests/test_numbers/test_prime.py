from cryptography318.number.prime import *
from cryptography318.number.crypto_functions import *
import pytest


def test_primes_lt():
    primes = primerange(11)
    assert primes == [2, 3, 5, 7, 11]
    assert len(primes) == prime_pi(11)

    primes = primes_gen(11)
    n = 2
    for p in primes:
        assert p == n
        n = next_prime(n)


def test_is_prime():
    n = 2047
    assert not isprime(n)

    n = randprime(pow(10, 3))
    assert confirm_prime(n)


def i_th_prime(i):
    j = 1
    prime = next_prime(1)
    while j < i:
        prime = next_prime(prime)
        j += 1
    return prime


def temp():
    lower = pow(10, 15)
    upper = pow(10, 16)

    a = randprime(lower, upper)
    b = randprime(lower, upper)
    print(f"a: {a}, b: {b}")

    x = a * b
    print(f"x : {x}; {len(str(x))}")


if __name__ == '__main__':
    temp()
