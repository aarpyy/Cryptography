from cryptography318.prime import randprime
from cryptography318.factor import qs
import random

from pathlib import Path


test_root = Path(__file__).parent.absolute()


def test_siqs():
    random.seed(10)
    lower = pow(10, 12)
    upper = lower * 10
    a = randprime(lower, upper)
    b = randprime(lower, upper)
    n = a * b
    c = qs(n, fp=test_root.parent.joinpath("cryptography318", "data", "primes.txt"))
    assert c == a or c == b


def test_siqs_trials():
    random.seed(10)
    lower = pow(10, 10)
    upper = lower * 10
    a = randprime(lower, upper)
    b = randprime(lower, upper)
    n = a * b

    # Use a low enough F that we won't find a factor
    F = 200
    assert qs(n, F=F, fp=test_root.parent.joinpath("cryptography318", "data", "primes.txt")) is None
