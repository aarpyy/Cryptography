import os
import sys
from random import seed

import pytest

from cryptography318.factor.elliptic import ecm_mont, gen_weierstrass_point, WeierstrassPoint, \
    ecm_mont_phase1, ecm_weierstrass, ecm_mont_basic, ecm
from cryptography318.prime import randprime


def test_weierstrass_point():
    seed(318)
    e = 7

    n = randprime(pow(10, e), pow(10, e + 1))

    P = gen_weierstrass_point(n)
    b = (pow(P.y, 2, n) - pow(P.x, 3, n) - P.a * P.x) % n

    Q = P + P

    # Assert Q is on curve
    assert (pow(Q.y, 2, n) - pow(Q.x, 3, n) - Q.a * Q.x) % n == b

    # Assert Q is not infinity
    assert Q is not WeierstrassPoint.infinity

    # Add infinity
    assert P + WeierstrassPoint.infinity == P

    # Multiply
    assert P * 2 == P + P

    # Multiply by zero
    assert P * 0 is WeierstrassPoint.infinity

    # Multiply by one
    assert P * 1 == P

    # Multiply by negative
    assert P * -1 == -P

    # Definition of identity
    assert P + -P is WeierstrassPoint.infinity


def test_ecm_weierstrass():
    seed(318)
    e = 5   # Weierstrass is slow so choose lower exponent

    a = randprime(pow(10, e), pow(10, e + 1))
    b = randprime(pow(10, e), pow(10, e + 1))

    n = a * b * b

    # Assert result divides n
    assert n % ecm_weierstrass(n) == 0

    # Assert no errors on verbose
    sys.stdout = open(os.devnull, "w")
    ecm_weierstrass(n, verbose=True)


def test_montgomery_point():
    seed(318)
    e = 7

    n = randprime(pow(10, e), pow(10, e + 1))

    P = ecm_mont_phase1(n, use_z1=False)

    Q = P.double()
    diff = P    # Difference between Q and P is P

    R = P.add(Q, diff)  # R = P + 2P
    assert R == P.ladder(3)

    assert P.double().double() == P.ladder(4)

    P = ecm_mont_phase1(n, use_z1=True)

    Q = P.double()
    diff = P    # Difference between Q and P is P

    R = P.add(Q, diff)  # R = P + 2P
    assert R == P.ladder(3)

    assert P.double().double() == P.ladder(4)


def test_ecm_mont_basic():
    seed(318)
    e = 5   # The basic algorithm is slow so choose lower exponent

    a = randprime(pow(10, e), pow(10, e + 1))
    b = randprime(pow(10, e), pow(10, e + 1))

    n = a * b * b

    # Assert result divides n
    assert n % ecm_mont_basic(n) == 0


def test_ecm_mont():
    seed(318)
    e = 10  # This should be much faster

    a = randprime(pow(10, e), pow(10, e + 1))
    b = randprime(pow(10, e), pow(10, e + 1))

    n = a * b * b

    # Assert result divides n
    assert n % ecm_mont(n, use_z1=False) == 0
    assert n % ecm_mont(n, use_z1=True) == 0

    # Assert no errors on verbose
    sys.stdout = open(os.devnull, "w")
    ecm_mont(n, use_z1=False, verbose=True)
    ecm_mont(n, use_z1=True, verbose=True)


def test_ecm():
    seed(318)
    e = 12

    a = randprime(pow(10, e), pow(10, e + 1))
    b = randprime(pow(10, e), pow(10, e + 1))

    n = a * b * b

    # Assert result divides n
    assert n % ecm(n, use_z1=False) == 0
    assert n % ecm(n, use_z1=True) == 0

    # Assert no errors on verbose
    sys.stdout = open(os.devnull, "w")
    ecm(n, use_z1=False, verbose=True)
    ecm(n, use_z1=True, verbose=True)


# @pytest.mark.skip(reason="Takes too long")
def test_specific_factor():
    n = 1362605012056345682341996237785766428743
    B1 = 10000
    for _ in range(10):
        value = ecm(n, B1=B1, verbose=True)
        if value is not None:
            print(value)
            return
        B1 *= 5

    assert False
