from cryptography318.linalg import (
    dot, flatten, matmul, rref, kernel, binary_kernel, det, eigvec, eigvals, minor, char_poly, matrix_equals
)
from cryptography318.linalg.linalg import matrix_copy, matrix_slice, identity_matrix, make_pivot, strmat
from cryptography318.linalg.matrix import Matrix
import pytest


def test_print_matrix():
    print(strmat([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))


def test_matrix_equals():
    assert matrix_equals([[1, 2, 3], [4, 5, 6]], ((1, 2, 3), (4, 5, 6)))
    assert matrix_equals([[1, 2, 3], [4, [5, 7], 6]], ((1, 2, 3), (4, [5, 7], 6)))
    assert not matrix_equals([[1, 2, 3], [4, 5, 6]], ((1, 2, 3), (4, 6)))
    assert not matrix_equals([[1, 2, 3], [4, 5, 6]], ((1, 2, 3), (4, 7, 6)))
    assert matrix_equals([], [])


def test_matrix_copy():
    a = [[1, 2, 3], [4, 5, 6]]
    b = matrix_copy(a)
    assert matrix_equals(b, a)
    b[0] = [2, 2, 3]
    assert a[0][0] == 1 and b[0][0] == 2


def test_matrix_slice():
    a = [[1, 2, 3], [4, 5, 6]]
    b, c = matrix_slice(a, 1)
    assert matrix_equals(b, [[1], [4]])
    assert matrix_equals(c, [[2, 3], [5, 6]])


def test_identity_matrix():
    I3 = identity_matrix(3)
    assert matrix_equals(I3, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    I0 = identity_matrix(0)
    assert matrix_equals(I0, [])


def test_make_pivot():
    a = [1, 2, 3, 4, 5]
    assert make_pivot(a, 1) == [e / 2 for e in a]
    assert make_pivot(a, 0) == a


def test_dot():
    a = [1, 2, 3, 4, 5]
    b = [1, 1, 1, 1, 1]
    assert dot(a, b) == 15
    b = [1, 1, 1, 1, -1]
    assert dot(a, b) == 5


def test_flatten():
    assert flatten([[1, 2, [[[3, 4, [[[10]]]]]]]]) == [1, 2, 3, 4, 10]


@pytest.mark.parametrize('args', [
    ([[1, 2, 3], [4, 5, 6]],
     [[1, 2], [3, 4], [5, 6]],
     [[22, 28], [49, 64]])
])
def test_matmul(args):
    a, b, c = args
    # Procedural style
    assert matrix_equals(matmul(a, b), c)
    # And class based
    b = Matrix(b)
    # Works with just __rmatmul__
    assert a @ b == c
    a = Matrix(a)
    # And __matmul__
    assert a @ b == c


def test_rref():
    from fractions import Fraction
    b = [[1, 0, -0.5], [0, 1, Fraction(-7) / 6]]
    a = [[Fraction(x) for x in row] for row in [[-2, 0, 1], [1, 3, -4]]]
    assert matrix_equals(rref(a), b)

    a = Matrix(a)
    assert a.rref() == b

    a = [[1, 2, 1, -1], [2, 4, 1, -1], [1, 2, 0, 0]]
    b = [[1, 2, 0, 0], [0, 0, 1, -1], [0, 0, 0, 0]]
    assert matrix_equals(rref(a), b)

    a = Matrix(a)
    assert a.rref() == b


def test_kernel():
    a = [[1, 2, 1, -1], [2, 4, 1, -1], [1, 2, 0, 0]]
    b = kernel(a)
    assert matrix_equals(b, [[-2, 1, 0, 0], [0, 0, 1, 1]])
    for x in b:
        for y in a:
            assert dot(x, y) == 0

    from fractions import Fraction
    a = [[Fraction(x) for x in row] for row in [[6.0, 2, -2], [2, 1.0, -2], [-2, -2, 6.0]]]
    assert kernel(a) == [[-1, 4, 1]]


def test_binary_kernel():
    a = [[1, 2, 1, -1], [2, 4, 1, -1], [1, 2, 0, 0]]
    b = binary_kernel(a)
    assert matrix_equals(b, [[0, 1, 0, 0], [0, 0, 1, 1]])
    for x in b:
        for y in a:
            assert dot(x, y) % 2 == 0


def test_det():
    a = [[3, -1, 4], [-1, 5, -9], [2, -6, 5]]
    assert det(a) == -90


def test_eigvec():
    from fractions import Fraction
    a = [[Fraction(x) for x in row] for row in [[0, 2, -2], [2, -5, -2], [-2, -2, 0]]]
    assert eigvec(a, values=list(int(x.real) if isinstance(x, complex) else int(x) for x in eigvals(a))) == {
        -6: [-1, 4, 1], -2: [1, 0, 1], 3: [-1, -0.5, 1]}


def test_eigvals():
    a = [[0, 2, -2], [2, -5, -2], [-2, -2, 0]]
    assert all(isinstance(x, float) for x in eigvals(a)) and set(eigvals(a)) == {-6.0, -2.0, 3.0}


def test_minor():
    a = [[3, -1, 4], [-1, 5, -9], [2, -6, 5]]
    assert minor(a) == -90


def test_char_poly():
    a = [[0, 2, -2], [2, -5, -2], [-2, -2, 0]]
    assert str(char_poly(a)) == "-x*(-x*(-x - 5) - 4) + 8*x + 36"


@pytest.mark.skip("All tests, only run specifically")
def test_all():
    test_matrix_slice()
    test_matrix_equals()
    test_matrix_copy()
    test_identity_matrix()
    test_make_pivot()
    test_dot()
    test_flatten()
    test_matmul()
    test_rref()
    test_kernel()
    test_binary_kernel()
    test_det()
    test_eigvals()
    test_eigvec()
    test_minor()
    test_char_poly()
