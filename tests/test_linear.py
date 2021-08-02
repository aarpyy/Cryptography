from cryptography318.linear_algebra_deprecated import (augRREF, IsSolvable, Solve, IsConsistent, MatrixEquals, RREF, IsRREF,
                                                       augIsRREF, _augRREF_old)
from cryptography318.matrix_deprecated import RandomMatrix, augmentMatrix, MultiplyMatrix, separateMatrix
from cryptography318.crypto_functions import QuadraticSieve, FactorInt, _factorPerfectSquare
from cryptography318.prime import RandomPrime, IsPrime
import textwrap
from cryptography318.linear_algebra import *
import numpy, random


def test_linear_solve(it=500):
    for _ in range(it):
        m = Matrix(rand=True, aug=True)
        m = m.rref()
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
    test_rref(50)
    test_linear_solve(50)


if __name__ == '__main__':
    testAllTests()
