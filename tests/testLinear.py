from cryptography318.linear_algebra import (augRREF, IsSolvable, Solve, IsConsistent, MatrixEquals, RREF, IsRREF,
                                            augIsRREF, _augRREF_old)
from cryptography318.matrix import RandomMatrix, augmentMatrix, MultiplyMatrix, separateMatrix
from cryptography318.crypto_functions import QuadraticSieve, FactorInt, _factorPerfectSquare
from cryptography318.prime import RandomPrime, IsPrime
import textwrap
import numpy, random


def testLinearSolve(it=500, count=False):
    acceptable, solved = 0, 0
    for _ in range(it):
        m1 = augRREF(RandomMatrix())
        if IsConsistent(m1) and IsSolvable(m1, aug=True):
            acceptable += 1
            res = Solve(m1)
            if res is not False:
                solved += 1
            if not count:
                assert res is not False

    if count:
        print(f"\nLinearSolve test tried to solve {it} random linear systems.")
        print(textwrap.indent(f"->Consistent and solvable systems: {acceptable}\n"
                              f"->Solved: {solved}\n", ' ' * 3))


def testRREF(it=500):
    for _ in range(it):
        m = RandomMatrix()
        m[0][0] = 0
        a = augRREF(m)
        assert augIsRREF(a)
        if IsConsistent(a) and IsSolvable(a, aug=True):
            assert type(Solve(a)) is numpy.ndarray


def testMatrixEquals():
    m2 = numpy.array([[44],
                      [25],
                      [26],
                      [7],
                      [45]])
    res = numpy.array([[43.99438],
                       [24.99189],
                       [25.98973],
                       [6.98902],
                       [45.00247]])
    assert MatrixEquals(m2, res)


def testAllTests():
    import time
    # testLinearSolve(it=1000, count=False)
    # testRREF()
    # testMatrixEquals()
    factors = QuadraticSieve(779, 7)
    power = 6
    p = RandomPrime(pow(10, power), pow(10, power + 1)) * RandomPrime(pow(10, power + 1), pow(10, power + 2))
    print(p)
    start2 = time.time()
    f = FactorInt(p)
    print(f"Factor Int took {time.time() - start2:.2f}s to factor {p + 1} into:\n{f}")
    for n in f:
        print(IsPrime(n))
    start1 = time.time()
    q = QuadraticSieve(p)
    print(f"Quadratic Sieve took {time.time() - start1:.2f}s to factor {p + 1} into:\n{q}")



if __name__ == '__main__':
    testAllTests()
