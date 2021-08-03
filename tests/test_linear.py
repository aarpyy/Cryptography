from cryptography318.linear_algebra_deprecated import (augRREF, IsSolvable, Solve, IsConsistent, MatrixEquals, RREF, IsRREF,
                                                       augIsRREF, _augRREF_old)
from cryptography318.matrix_deprecated import RandomMatrix, augmentMatrix, MultiplyMatrix, separateMatrix
from cryptography318.crypto_functions import QuadraticSieve, FactorInt, _factorPerfectSquare
from cryptography318.prime import RandomPrime, IsPrime, NextPrime
import textwrap
from cryptography318.linear_algebra import *
import numpy, random


def test_linear_solve(it=500):
    for _ in range(it):
        m = Matrix(rand=True, aug=True)
        print(m)
        m = m.rref()
        print(m)
        if m.is_solvable():
            assert m.solve()


def test_rref(it=500):
    for _ in range(it):
        m = RandomMatrix()
        m[0][0] = 0
        a = augRREF(m)
        assert augIsRREF(a)
        if IsConsistent(a) and IsSolvable(a, aug=True):
            assert type(Solve(a)) is numpy.ndarray


def testAllTests():
    # test_rref(50)
    test_linear_solve(1)


if __name__ == '__main__':
    # testAllTests()
    m = Matrix(rows=4, cols=4, identity=True)
    x = m.array
    x[0][1] = 5
    print(x)
    print(x == 1)
    # m[1][0] = 1
    # m[2][3] = 1
    # m[3][3] = 0
    # print(m)
    # m.is_rref1()
    # x = where(m == 1)
    # print(x)
    # for e in zip(x[0], x[1]):
    #     print(e)

    # print(m == 1)
    # print(type(m[0][0]))
    # m /= 5
    # print(m)
