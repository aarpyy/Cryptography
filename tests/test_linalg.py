from cryptography318.tools import *
from cryptography318.linalg import *
import pytest
import numpy


def test_rref():
    for _ in range(500):
        a = Matrix(rand=True)
        a = a.rref()
        assert a.is_rref() and a.is_rref_old()


if __name__ == '__main__':
    test_rref()
    a = Matrix(rand=True)
    b = a.copy()
    print(b == a)
