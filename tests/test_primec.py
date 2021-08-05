from pathlib import Path
from ctypes import CDLL
from cryptography318.prime import IsPrime
from random import randrange
import pytest


@pytest.mark.skip
def test_cprime():
    libname = Path().absolute() / "cryptography318/libprime.so"
    c_lib = CDLL(libname)
    for _ in range(50):
        n = randrange(3, pow(2, 30))
        # print(n, bool(c_lib.is_prime(n)))



test_cprime()
