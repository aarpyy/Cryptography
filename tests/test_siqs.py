from cryptography318.prime import randprime
from cryptography318.siqs import siqs

from pathlib import Path


test_root = Path(__file__).parent.absolute()


def test_siqs_34():
    a = randprime(pow(10, 17), pow(10, 18))
    b = randprime(pow(10, 17), pow(10, 18))
    c = siqs(a * b, fp=str(test_root.parent.joinpath("cryptography318").joinpath("primes.txt")), loud=False)
    assert c == a or c == b
