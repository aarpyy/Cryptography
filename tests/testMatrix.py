from cryptography318.matrix_deprecated import MakeMatrix, MultiplyMatrix, Transpose
from cryptography318.linear_algebra_deprecated import InvertMatrix, ChangeBasisMap, MatrixEquals
from cryptography318.linear_algebra import *
import numpy, random
import pytest

S = numpy.array([[2, -1],
                 [5, -3],
                 [-3, 2]])

C = numpy.array([[1, 1, 0],
                 [1, 0, 1],
                 [0, 1, 1]])

B = numpy.array([[2, 3],
                 [3, 5]])


def testChangeBasis():
    test1 = ChangeBasisMap(S, in_basis=numpy.array([[1, 0], [0, 1]]),
                           out_basis=numpy.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
    test2 = ChangeBasisMap(S, in_basis=numpy.array([[1, 0], [0, 1]]))
    test3 = ChangeBasisMap(S, out_basis=numpy.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
    S1 = LinearMap(S)
    C1 = Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    B1 = Matrix([[1, 0], [0, 1]])
    test4 = S1.change_basis(in_basis=B1, out_basis=C1)
    test5 = S1.change_basis(in_basis=B1)
    test6 = S1.change_basis(out_basis=C1)

    m = Matrix(rows=5, cols=5, identity=True, aug=True)
    m[4][4] = 0
    m[3][4] = 3
    m[2][2] = 2
    m[1][3] = 5
    m[4][4] = 1
    m[4][3] = 1

    test_map = LinearMap([[3, 2],
                          [-1, 0]])
    vector = Matrix([[1],
                     [1]])
    value = 2
    mat = Matrix([[3, 2],
                  [-1, 0]])


def testRankNull():
    for _ in range(5000):
        x = random.randrange(2, 10)
        m = Matrix(rows=x, cols=x, rand=True)
        assert m.rank() + m.null() == m.dimension()


def testRREF():
    for _ in range(5000):
        x = random.randrange(2, 10)
        m = Matrix(rows=x, cols=x, rand=True)
        I_x = Matrix(rows=x, cols=x, identity=True)
        assert m.rref() == I_x


def testAllTests():
    testChangeBasis()
    testRankNull()
    testRREF()


def homework_9():
    """Problem 3.5.8"""
    b = Matrix([[1, 1],
                [1, -1]])
    c = Matrix([[1, 1, 0],
                [1, 0, 1],
                [0, 1, 1]])
    a1 = Matrix([[-2],
                 [1]])
    a2 = Matrix([[2],
                 [3]])
    a3 = Matrix([[1],
                 [-2],
                 [3]])
    a4 = Matrix([[3],
                 [2],
                 [1]])
    T = LinearMap([[0, 1, 0],
                [1, 0, 0]])

    print("a)\ni)\n1/2 *")
    print(2 * a1.change_basis(b))
    print("\nii)\n1/2 *")
    print(2 * a2.change_basis(b))
    print("\niii)")
    print(a3.change_basis(c))
    print("\niv)")
    print(a4.change_basis(c))
    print("\nb)\nv)")
    print(T.change_basis(in_basis=c))
    print("\nvi)\n1/2 *")
    print((2 * T.change_basis(out_basis=b)))
    print("\nvii)\n1/2 *")
    print((2 * T.change_basis(in_basis=c, out_basis=b)))


if __name__ == '__main__':
    # testAllTests()
    homework_9()
