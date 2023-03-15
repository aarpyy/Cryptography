from cryptography318.prime import randprime
from cryptography318.factor.siqs import siqs
import random

from pathlib import Path


test_root = Path(__file__).parent.absolute()


def test_siqs_34():
    random.seed(10)
    lower = pow(10, 12)
    upper = lower * 10
    a = randprime(lower, upper)
    b = randprime(lower, upper)
    n = a * b
    print(f"Factors of {n}: {{a: {a}, b: {b}}}")
    c = siqs(n, fp=str(test_root.parent.joinpath("cryptography318").joinpath("prime/primes.txt")), loud=True)
    assert c == a or c == b
