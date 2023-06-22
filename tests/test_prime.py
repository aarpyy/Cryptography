from cryptography318.prime import sqrt_mod, randprime, quadratic_residue, lift_sqrt, next_prime, isprime, prev_prime, get_details
from cryptography318.prime.prime import primesieve
from random import randrange, seed
import json


def test_next_prime():
    seed(318)
    l = pow(10, 4)
    for _ in range(50):
        n = randrange(l, l * 10)
        p = next_prime(n)
        assert isprime(p) and p >= n

        # Everything in between is not prime
        for i in range(n + 1, p):
            assert not isprime(i), f"{n} <= {i} <= {p} is prime"


def test_prev_prime():
    seed(318)
    l = pow(10, 4)
    for _ in range(50):
        n = randrange(l, l * 10)
        p = prev_prime(n)
        assert isprime(p) and p <= n

        # Everything in between is not prime
        for i in range(p + 1, n, -1):
            assert not isprime(i), f"{n} <= {i} <= {p} is prime"


def test_sqrt_mod():
    seed(318)

    for _ in range(50):
        p = randprime(pow(10, 3), pow(10, 4))
        a = randrange(2, p)

        # Need a value that is a quadratic residue
        while not quadratic_residue(a, p):
            a = randrange(2, p)

        s = sqrt_mod(a, p)
        assert s is not None and pow(s, 2, p) == a


def test_lift_sqrt():
    seed(318)
    e = pow(10, 3)
    for _ in range(50):
        p = randprime(e, e * 10)

        # Get random square root
        root = randrange(p // 2, p)
        a = pow(root, 2, p)
        for i in range(1, 4):
            m = p ** i
            root = lift_sqrt(root, a, m, p)
            assert pow(root, 2, m) == a % m


def test_primesieve():
    seed(318)
    n = 200
    k = 5
    primesieve.extend(n)

    for _ in range(200):
        value = randrange(2, n)
        result = primesieve.search(value)
        if isprime(value):
            assert isinstance(result, int) and primesieve[result] == value
        else:
            assert isinstance(result, tuple) and len(result) == 2 and all(isinstance(r, int) for r in result) and \
                   primesieve[result[0]] < value < primesieve[result[1]], f"malformed result: {result} for {value}"

        # Make sure to remove duplicates
        values = sorted([*set([randrange(2, n) for _ in range(k)])])
        results = primesieve.search(*values)

        assert len(values) == len(results)

        for value, result in zip(values, results):
            if isprime(value):
                assert isinstance(result, int) and primesieve[result] == value
            else:
                assert isinstance(result, tuple) and len(result) == 2 and all(isinstance(r, int) for r in result) and \
                       primesieve[result[0]] < value < primesieve[result[1]], f"malformed result: {result} for {value}"


def test_get_details():
    seed(318)
    for _ in range(50):
        e = randrange(3, 10)
        n = randrange(pow(10, e), pow(10, e + 1)) ^ 1
        isprime(n, details=True)

        details = get_details()

        # Assert at least one description is provided
        assert any(v is not None for v in details.values())
