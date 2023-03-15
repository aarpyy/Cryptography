from cryptography318.prime import sqrt_mod, randprime, quadratic_residue, lift_sqrt
from cryptography318.prime.prime import primesieve
from random import randrange


def test_sqrt_mod():
    for _ in range(50):
        p = randprime(pow(10, 3), pow(10, 4))
        a = randrange(2, p)
        while not quadratic_residue(a, p):
            a = randrange(2, p)

        s = sqrt_mod(a, p)
        assert s is not None and pow(s, 2, p) == a


def test_lift_sqrt():
    for _ in range(50):
        p = randprime(pow(10, 3), pow(10, 4))
        root = randrange(p // 2, p)
        a = pow(root, 2, p)
        assert pow(root, 2, p) == a % p
        m = p
        for __ in range(3):
            root = lift_sqrt(root, a, m, p)
            m *= p
            assert pow(root, 2, m) == a % m


def test_primesieve():
    primesieve.extend(200)
    print(primesieve)
    r = primesieve.search(4, 53, 59, 62, 71)
    print(r)
    print(primesieve[r[3][0]])
    print(primesieve[r[3][1]])


def test_temp():
    from math import lcm, gcd
    l = pow(2, 40)
    u = l << 2
    p = randprime(l, u)
    q = randprime(l, u)
    n = p * q
    y = lcm(p - 1, q - 1)
    e = pow(2, 16) + 1
    while gcd(e, y) != 1:
        e += 1
    d = pow(e, -1, y)

    m = randrange(pow(2, 25))

    c = pow(m, d, n)
    assert pow(c, e, n) == m % n

    print(f"n: {n}; p: {p}; q: {q}; e: {e}; d: {d}")
