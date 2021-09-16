from cryptography318.linalg.linalg import Matrix, is_binary_matrix
import numpy
import pytest


def is_float(n):
    return isinstance(n, (float, numpy.float16, numpy.float32, numpy.float64))


def is_int(n):
    return isinstance(n, (int, numpy.int16, numpy.int32, numpy.int64))


def test_constructor():
    m = Matrix([[1, 2, 3],
                [3, 2, 1]])

    assert isinstance(m, Matrix)
    assert isinstance(m.array, numpy.ndarray)
    assert m.augmented is False

    # construction with np.array
    np_array = numpy.array([[1, 2, 3],
                            [3, 2, 1]])
    m = Matrix(np_array)
    assert isinstance(m, Matrix)
    assert isinstance(m.array, numpy.ndarray)
    assert m.augmented is False

    # construction of identity matrix
    m = Matrix(rows=3, cols=3, identity=True)
    assert isinstance(m, Matrix)
    assert isinstance(m.array, numpy.ndarray)
    assert m == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

    # construction of identity matrix with diag != 1
    m = Matrix(rows=3, cols=3, identity=5)
    assert isinstance(m, Matrix)
    assert isinstance(m.array, numpy.ndarray)
    assert m == [[5, 0, 0], [0, 5, 0], [0, 0, 5]]

    # construction of random matrix
    m = Matrix(rand=True)
    assert isinstance(m, Matrix)
    assert isinstance(m.array, numpy.ndarray)

    # construction of null matrix
    m = Matrix(rows=3, cols=3)
    assert isinstance(m, Matrix)
    assert isinstance(m.array, numpy.ndarray)
    assert m == [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

    # construction with non-list-type object as array
    with pytest.raises(TypeError) as exc_info:
        Matrix(5)
    assert "type numpy.ndarray or list" in str(exc_info.value)

    # null construction
    with pytest.raises(ValueError):
        Matrix()

    # construction of row vector returns nested row vector
    m = Matrix([1, 2, 3])
    assert isinstance(m[0], numpy.ndarray)


def test_make_pivot():
    m = Matrix(rows=4, cols=4, identity=3)
    m[1][0] = 6
    m.make_pivot(1)
    assert m[1][0] == 1

    m = Matrix(rows=4, cols=4, identity=3)
    m[1][0] = 6
    m.make_pivot(1, col=1)
    assert m[1][0] == 2 and m[1][1] == 1


def test_operators():
    m = Matrix([[1, 2, 3],
                [3, 2, 1]])
    a = 5 * m
    assert a[0][0] == 5

    # can matrix convert to float when using float division
    a = m / 2
    assert is_float(a[0][0])

    a = m // 2
    assert is_int(a[0][0])

    a = m + 2
    assert a[0][0] == 3

    a = m - 2
    assert a[0][0] == -1

    a = m % 2
    assert is_binary_matrix(a)

    I = Matrix(rows=3, cols=3, identity=True)
    inv = pow(I, -1)
    assert inv == I

    m = Matrix(rows=3, cols=3, rand=True)
    a = m * m * m
    assert pow(m, 3) == a

    with pytest.raises(TypeError) as exc_info:
        a = 5 - m
    assert "subtraction unsupported" in str(exc_info.value)


def test_reset_type():
    b = numpy.array([[2.0, 3.5, 4.6],
                     [5.4, 3, 9]])
    c = numpy.array([[2, 3., 4.],
                     [5., 3, 9]])

    assert is_float(c[0][0]) and is_float(b[0][0])

    # reset type is called when input of matrix is list-type
    m = Matrix(c)
    assert is_int(m[0][0])

    # some values cannot be converted to ints in b, so should be float
    m = Matrix(b)
    assert is_int(m[0][0]) and is_float(m[0][1])


def test_mod():
    b = Matrix([[-32, 35, 10, -12],
                [-5, -26, 26, 11],
                [-37, -48, -20, 28],
                [-46, -26, 40, 25],
                [-43, -38, -33, -14],
                [-22, 16, -33, 24]], aug=True, mod=13)

    b.reduce_mod()

    # after reducing everything should be positive integer within 0 and the row modulus
    for i in range(len(b)):
        for j, n in enumerate(b[i]):
            assert 0 <= n < b.mod and isinstance(n, (int, numpy.int64))

    possible_pivots = b.find_invertible()

    # the third row in transpose is the final column in standard form, which is augmented list of solutions
    # and there should be no possible pivot there
    assert len(numpy.where(possible_pivots.transpose()[3] == 1)[0]) == 0


if __name__ == '__main__':
    test_operators()

