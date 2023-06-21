import json
from pathlib import Path
from random import seed, randrange

import pytest

from cryptography318.factor import siqs, factor, pollard_p1, pollard_rho_factor, get_details
from cryptography318.prime import randprime

test_root = Path(__file__).parent.absolute()


def test_siqs():
    a = randprime(pow(10, 17), pow(10, 18))
    b = randprime(pow(10, 17), pow(10, 18))
    file = test_root.parent.joinpath("cryptography318").joinpath("primes.txt")
    if file.is_file():
        c = siqs(a * b, fp=str(file), loud=False)
    else:
        c = siqs(a * b, loud=False)
    assert c == a or c == b


def test_rho():
    a = randprime(pow(10, 12), pow(10, 13))
    b = randprime(pow(10, 12), pow(10, 13))
    c = pollard_rho_factor(a * b)
    assert c == a or c == b


def test_p1():
    a = randprime(pow(10, 12), pow(10, 13))
    b = randprime(pow(10, 12), pow(10, 13))
    c = randprime(pow(10, 2), pow(10, 3))
    d = pollard_p1(a * b * c)
    assert d in (a, b, c)


def test_factor():
    seed(318)
    e = 10
    for _ in range(10):
        a = randprime(pow(10, e), pow(10, e + 1))
        b = randprime(pow(10, e), pow(10, e + 1))
        c = randprime(pow(10, 2), pow(10, 3))
        factors = factor(a * b * b * c * c)
        assert factors == {a: 1, b: 2, c: 2}


def test_get_details():
    seed(318)
    e = 10
    for _ in range(10):
        n = randrange(pow(10, e), pow(10, e + 1))
        factor(n, details=True)
        details = get_details()
        details_obj = json.loads(details)

        # At least one of the values must be literal True
        assert any(v is True for v in details_obj.values())

        a = randprime(pow(10, e), pow(10, e + 1))
        b = randprime(pow(10, e), pow(10, e + 1))
        n = a * b * b
        factor(n, details=True)
        details = get_details()
        details_obj = json.loads(details)

        # At least one of the values must be literal True
        assert any(v is True for v in details_obj.values())


@pytest.mark.skip()
def test_all():
    test_siqs()
    test_rho()
    test_p1()
    test_factor()


if __name__ == "__main__":
    test_all()
