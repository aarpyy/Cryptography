from cryptography318.linear_algebra_deprecated import (augRREF, IsSolvable, Solve, IsConsistent, MatrixEquals, RREF, IsRREF,
                                                       augIsRREF, _augRREF_old)
from cryptography318.matrix_deprecated import RandomMatrix, augmentMatrix, MultiplyMatrix, separateMatrix
from cryptography318.crypto_functions import QuadraticSieve, FactorInt, _factorPerfectSquare
from cryptography318.prime import RandomPrime, IsPrime, NextPrime
import textwrap
from cryptography318.linear_algebra import *
import numpy, random, pytest


@pytest.mark.skip
def test_change_basis():
    S = LinearMap([[1, 0],
                   [0, -1]])
    B = Matrix([[1, -2],
                [2, 1]]).invert()
    T = Matrix([[-.6, .8],
                [.8, .6]])
    assert S.in_basis(B) == T


@pytest.mark.skip
def test_ranknull(it=500):
    for _ in range(it):
        x = random.randrange(2, 10)
        m = Matrix(rows=x, cols=x, rand=True)
        assert m.rank() + m.null() == m.dimension()


@pytest.mark.skip
def test_rref(it=500):
    for _ in range(it):
        x = random.randrange(2, 10)
        m = Matrix(rows=x, cols=x, rand=True)
        assert m.rref().is_rref()


@pytest.mark.skip
def test_basis():
    A = LinearMap([[2, 2],
                   [1, 3]])
    print(A.is_eigen_value(1))


@pytest.mark.skip
def test_linear_solve(it=500):
    for _ in range(it):
        m = Matrix(rand=True, aug=True)
        print(m)
        m = m.rref()
        print(m)
        if m.is_solvable():
            assert m.solve()


@pytest.mark.skip
def test_rref(it=500):
    for _ in range(it):
        m = RandomMatrix()
        m[0][0] = 0
        a = augRREF(m)
        assert augIsRREF(a)
        if IsConsistent(a) and IsSolvable(a, aug=True):
            assert type(Solve(a)) is numpy.ndarray


@pytest.mark.skip
def testAllTests():
    test_rref(50)
    test_linear_solve(1)


if __name__ == '__main__':
    testAllTests()
