from cryptography318.crypto_functions import *
from cryptography318.prime import *


def test_bsgs():
    for _ in range(50):
        g = randrange(3, 10)
        p = RandomPrime(pow(10, 2), pow(10, 3))
        x = randrange(2, p - 1)
        h = pow(g, x, p)

        assert pow(g, baby_step_giant_step(g, h, p), p) == h
