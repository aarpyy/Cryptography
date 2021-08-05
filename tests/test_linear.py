from cryptography318.linear_algebra_deprecated import (augRREF, IsSolvable, Solve, IsConsistent, MatrixEquals, RREF, IsRREF,
                                                       augIsRREF, _augRREF_old)
from cryptography318.matrix_deprecated import RandomMatrix, augmentMatrix, MultiplyMatrix, separateMatrix
from cryptography318.crypto_functions import QuadraticSieve, FactorInt, _factorPerfectSquare
from cryptography318.prime import RandomPrime, IsPrime, NextPrime
import textwrap
from cryptography318.linear_algebra import *
import numpy, random, pytest


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
    # test_rref(50)
    test_linear_solve(1)


if __name__ == '__main__':
    # testAllTests()
    mat = Matrix(rows=3, cols=2)
    mat[1][0] = 1
    arr = numpy.array([[1, 0, 2]]).astype(numpy.float64)

    m = array_mod([15, 6, 9, 15], mod=18)
    a = numpy.array([[1, 5, 7, 3, 4]])
    x = [1, 2, 3, 5, 3]
    c = [1, 2, 3]
    d = [1, 0, 2]
    print(m, m.mod)
    m = m - 10
    print(m, m.mod)


