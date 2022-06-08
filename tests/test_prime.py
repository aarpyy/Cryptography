from cryptography318.prime import sqrt_mod, randprime, quadratic_residue, lift_sqrt
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
