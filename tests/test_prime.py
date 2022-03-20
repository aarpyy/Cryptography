from cryptography318.prime import sqrt_mod, randprime, quadratic_residue
from random import randrange


def test_sqrt_mod():
    for _ in range(50):
        p = randprime(pow(10, 3), pow(10, 4))
        a = randrange(2, p)
        while not quadratic_residue(a, p):
            a = randrange(2, p)

        s = sqrt_mod(a, p)
        assert s is not None and pow(s, 2, p) == a
