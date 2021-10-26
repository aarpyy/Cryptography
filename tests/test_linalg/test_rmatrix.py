from cryptography318.linalg.matrix import *
from cryptography318.linalg.arrayabc import *
from timeit import timeit
import numpy as np
from sympy import Symbol
import sympy as sy
from typing import *
from fractions import Fraction
from numbers import *
from cryptography318.core.tools import *
import operator


def test_rref():
    a = rmatrix([[1, 0, 0, 0, 3],
                 [0, 0, 1, 0, 0],
                 [0, 1, 0, 1, 1]])
    assert not a.is_rref

    a = rmatrix([[1, 0, 0, 0, 3],
                 [0, 0, 1, 0, 0],
                 [0, 0, 0, 1, 1]])
    assert a.is_rref

    a = rmatrix([[1, 0, 0, 0, 3],
                 [0, 0, 1, 0, 0],
                 [0, 0, 0, 0, 1]])
    assert not a.is_rref

    a = rmatrix([[1, 0, 0, 0, 3],
                 [0, 0, 1, 0, 0],
                 [0, 0, 0, 0, 1]], aug=True)
    assert a.is_rref


def test_consistent():
    a = rmatrix([[1, 0, 0, 0, 3],
                 [0, 0, 1, 0, 0],
                 [0, 0, 0, 0, 0]], aug=True)
    print(a.rank)


def test_transpose():
    _iters = pow(10, 5)
    a = rmatrix([[1, 2, 3], [4, 5, 6]])
    as_list = a.aslist
    nparr = np.array(as_list)
    time_py = timeit(lambda: list(map(list, zip(*as_list))), number=_iters)
    time_rm = timeit(lambda: a.transpose(), number=_iters)
    time_np = timeit(lambda: nparr.transpose(), number=_iters)
    print(f"time rmatrix: {time_rm:.3f}")
    print(f"time numpy: {time_np:.3f}")
    print(f"time python: {time_py:.3f}")


def test_operations():
    a = rmatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    b = [[2], [4], [6]]
    _iters = pow(10, 3)
    print(timeit(lambda: a * b, number=_iters))
    c = np.array([[Fraction(1), Fraction(2), Fraction(3)], [Fraction(4), Fraction(5), Fraction(6)]])
    print(timeit(lambda: c @ b, number=_iters))
    c = c.astype(int)
    print(type(c[0][0]))
    print(isinstance(c[0][0], Integral))


if __name__ == '__main__':
    pass
