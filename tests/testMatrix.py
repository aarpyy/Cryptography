from cryptography318.matrix import MakeMatrix, MultiplyMatrix, Transpose
from cryptography318.linear_algebra import InvertMatrix, ChangeBasisMap, MatrixEquals
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
    assert MatrixEquals(ChangeBasisMap(S, in_basis=numpy.array([[1, 0], [0, 1]]),
                                       out_basis=numpy.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])), S)
    assert MatrixEquals(ChangeBasisMap(S, in_basis=numpy.array([[1, 0], [0, 1]])), S)
    assert MatrixEquals(ChangeBasisMap(S, out_basis=numpy.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])), S)


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
