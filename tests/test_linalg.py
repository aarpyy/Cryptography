from cryptography318.tools import *
from cryptography318.linalg import *
from cryptography318.linear_algebra import Matrix as Matrix2
import pytest
import numpy
from cryptography318.fraction import *


def test_rref():
    for _ in range(500):
        a = Matrix(rand=True)
        a = a.rref()
        assert a.is_rref() and a.is_rref_old()


def test_remove_null():
    while False:
        a = Matrix(rand=True)
        if len(a) < 2 or len(a[0]) < 2 or a.remove_null(copy=True) != a:
            continue
        a[len(a) // 2] = Array([0] * len(a[0]))  # make null row
        b = a.transpose()
        b[len(b) // 2] = Array([0] * len(b[0]))  # make null column
        a = b.transpose()
        try:

            assert len(k := a.remove_null_column(copy=True)) == len(a)
            assert len(k[0]) == len(a[0]) - 1

            assert len(k := a.remove_null_row(copy=True)) == len(a) - 1
            assert len(k[0]) == len(a[0])

            assert len(k := a.remove_null(copy=True)) == len(a) - 1
            assert len(k[0]) == len(a[0]) - 1

            assert (c := a.copy()).mod == a.mod and c.augmented == a.augmented

        except:
            print(a)
            break

    a = Matrix([[0, 0],
                [0, 0]])
    assert not a.remove_null_row(copy=True)


def test_char_poly():
    a = Matrix([[-1, -3, 1],
                [3, 3, 1],
                [3, 0, 4]])
    x = Symbol('x')
    assert a.char_poly(sym=x) == -6 * x + (3 - x) * (4 - x) * (-x - 1) + 18


def test_eigvec():
    b = matrix([[1, 2, 3],
                [4, 5, 6],
                [6, 7, 8]])
    c = [[2],
         [3],
         [4]]
    print(e := b.eigvals())
    print(eig := b.eigvec(e))
    print(b * eig)


if __name__ == '__main__':
    # test_rref()
    # test_remove_null()
    # test_char_poly()
    a = matrix([[1, 0, -3, 0, 2, -8],
                [0, 1, 5, 0, -1, 4],
                [0, 0, 0, 1, 7, -9],
                [0, 0, 0, 0, 0, 0]])

    m = matrix([[-14, 2, 3], [4, -10, 6], [6, 7, -7]])
    # print(r := m.orthonormalize(steps=True))
    print(factor(72))

