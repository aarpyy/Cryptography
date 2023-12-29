from pathlib import Path
from random import seed, randrange
from time import time

import jsonschema
import pytest

from cryptography318 import ecm
from cryptography318.factor import qs, factor, pollard_pm1, pollard_rho_factor
from cryptography318.prime import randprime

test_root = Path(__file__).parent.absolute()


def test_compare_ecm():
    from sympy.ntheory.ecm import _ecm_one_factor

    seed(318)
    e = 7

    a = randprime(pow(10, e), pow(10, e + 1))
    b = randprime(pow(10, e), pow(10, e + 1))

    n = a * b * b

    B1 = 10000
    B2 = 100 * B1
    num_curves = 50

    start = time()
    factor_sympy = _ecm_one_factor(n, B1, B2, num_curves)
    time_sympy = time() - start

    print(f"sympy: {time_sympy}")

    start = time()
    factor_cr = ecm(n, use_z1=True, use_weierstrass=True)
    time_cr = time() - start

    print(f"cryptography318: {time_cr}")

    print(f"sympy: {factor_sympy}")
    print(f"cryptography318: {factor_cr}")


def test_compare_sympy():
    from sympy.ntheory import factorint

    seed(318)
    e = 11
    time_sympy = 0
    time_cr = 0

    a = randprime(pow(10, e), pow(10, e + 1))
    b = randprime(pow(10, e), pow(10, e + 1))
    n = a * b * b
    start = time()
    factors_sympy = factorint(n, verbose=True, use_rho=False, use_trial=False, use_pm1=False)
    time_sympy += time() - start

    start = time()
    details = {}
    factors_cr = factor(n, details=details, use_rho=False, use_pm1=False, use_siqs=False, verbose=True)
    time_cr += time() - start

    assert factors_cr == factors_sympy

    print(f"sympy: {time_sympy}")
    print(f"cryptography318: {time_cr}")
    print(details)


def test_siqs():
    a = randprime(pow(10, 17), pow(10, 18))
    b = randprime(pow(10, 17), pow(10, 18))
    file = test_root.parent.joinpath("cryptography318", "data", "primes.txt")
    if file.is_file():
        c = qs(a * b, fp=str(file), loud=False)
    else:
        c = qs(a * b, loud=False)
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
    d = pollard_pm1(a * b * c)
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


def test_details():
    schema = {
        "type": "object",
        "properties": {
            "methods": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "function": {"type": "string"},
                        "name": {"type": "string"},
                        "value": {},
                    },
                },
            },
            "error": {
                "type": "string",
            },
        },
        "required": ["methods"]
    }

    seed(318)
    e = 10
    n = randrange(pow(10, e), pow(10, e + 1))
    details = {}
    factor(n, details=details)

    jsonschema.validate(details, schema)
    assert "error" not in details

    a = randprime(pow(10, e), pow(10, e + 1))
    b = randprime(pow(10, e), pow(10, e + 1))
    n = a * b * b
    details = {}
    factor(n, details=details)

    jsonschema.validate(details, schema)
    assert "error" not in details

    details = {}
    factor(1.2, details=details)
    jsonschema.validate(details, schema)
    assert "error" in details


@pytest.mark.skip()
def test_all():
    test_siqs()
    test_rho()
    test_p1()
    test_factor()


if __name__ == "__main__":
    test_details()
