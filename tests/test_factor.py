import json
from pathlib import Path
from random import seed, randrange

import jsonschema
import pytest

from cryptography318.factor import siqs, factor, pollard_p1, pollard_rho_factor
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
                        "value": {
                            "description": "The value returned by the function",
                        },
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
