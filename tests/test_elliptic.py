from cryptography318.factor.elliptic import ecm_weierstrass, ecm_mont
from cryptography318.prime import randprime
import random


def test_ecm_mont():
    random.seed(10)
    failed = 0
    n = 10
    for _ in range(n):
        a = randprime(pow(10, 25), pow(10, 26))
        b = randprime(pow(10, 25), pow(10, 26))
        c = randprime(pow(10, 4), pow(10, 5))
        d = ecm_mont(a * b * c, B1=960, B2=5700)
        if d is None:
            failed += 1
            continue
        assert d in (a, b, c)
    assert failed < n, "Failed every time"


def test_ecm_weierstrass():
    random.seed(10)
    failed = 0
    n = 10
    for _ in range(n):
        a = randprime(pow(10, 7), pow(10, 8))
        b = randprime(pow(10, 7), pow(10, 8))
        f = ecm_weierstrass(a * b)
        if f is None:
            failed += 1
            continue
        assert f in (a, b)
    assert failed < n, "Failed every time"
