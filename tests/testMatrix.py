from cryptography318.matrix_deprecated import MakeMatrix, MultiplyMatrix, Transpose
from cryptography318.linear_algebra_deprecated import InvertMatrix, ChangeBasisMap, MatrixEquals
from cryptography318.linear_algebra import *
import numpy
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
    print(S1.map(C1))


def testAllTests():
    testChangeBasis()


"""Problem 5: a), b), c)"""

# a = ChangeBasisMap(S, in_basis=B)
# b = ChangeBasisMap(S, out_basis=C)
# c = ChangeBasisMap(S, in_basis=B, out_basis=C)
#
# print(f"[S]B->E\n{a}\n")
# print(f"[S]E->C\n{b}\n")
# print(f"[S]B->C\n{c}\n")

"""Problem 4: b), c)"""

# basis = numpy.array([[1, 1, 0], [-1, 2, 2], [1, -1, -1]])
# for _ in range(2):
#     m = MakeMatrix(3, 1)
#     print(MultiplyMatrix(InvertMatrix(basis), m))

if __name__ == '__main__':
    testAllTests()
