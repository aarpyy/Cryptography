from cryptography318.linear_algebra import *
from random import randrange
import pytest


def test_change_basis():
    S = LinearMap([[1, 0],
                   [0, -1]])
    B = Matrix([[1, -2],
                [2, 1]]).invert()
    T = Matrix([[-.6, .8],
                [.8, .6]])
    assert S.in_basis(B) == T


def test_ranknull(it=50):
    for _ in range(it):
        x = randrange(2, 10)
        m = Matrix(rows=x, cols=x, rand=True)
        assert m.rank() + m.null() == m.dimension()


def test_rref(it=50):
    for _ in range(it):
        x = randrange(2, 10)
        m = Matrix(rows=x, cols=x, rand=True)
        assert m.rref().is_rref()


def test_basis():
    A = LinearMap([[2, 2],
                   [1, 3]])
    assert A.is_eigen_value(1)


def test_linear_solve(it=50):
    for _ in range(it):
        m = Matrix(rand=True, aug=True)
        m = m.rref()
        if m.is_solvable():
            assert m.solve()


if __name__ == '__main__':
    test_change_basis()
    test_ranknull()
    test_rref()
    test_basis()
    test_linear_solve()


