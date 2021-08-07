from cryptography318.linear_algebra import *
import numpy
import pytest


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

    # construction of row vector returns layered row vector
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


def test_reset_type():
    b = numpy.array([[2.0, 3.5, 4.6],
                     [5.4, 3, 9]])
    c = numpy.array([[2, 3., 4.],
                     [5., 3, 9]])
    m = Matrix(b)
    a = Matrix(c)
    m.reset_type()
    a.reset_type()
    print(m)
    print(type(m[0][0]))
    print(a)
    print(type(a[0][0]))


@pytest.mark.skip
def test_mod():
    mods = [20] * 6
    mods[0] = 13
    b = Matrix([[-32, 35, 10, -12],
                [-5, -26, 26, 11],
                [-37, -48, -20, 28],
                [-46, -26, 40, 25],
                [-43, -38, -33, -14],
                [-22, 16, -33, 24]], mod=mods, aug=True)

    b.reduce_mod()
    print(b)
    # after reducing everything should be positive integer within 0 and the row modulus
    for i in range(len(b)):
        for j, n in enumerate(b[i]):
            assert 0 <= n < b.mod[i] and isinstance(n, (int, numpy.int64))

    possible_pivots = b.find_pivots_mod()
    print(possible_pivots)

    # the third row in transpose is the final column in standard form, which is augmented list of solutions
    # and there should be no possible pivot there
    assert len(numpy.where(possible_pivots.transpose()[3] == 1)[0]) == 0

    pivots = b.choose_pivots_mod()

    print(pivots)


def test_mod2():
    mods = [20] * 6
    mods[0] = 13
    m = Matrix([[6, 8, 9, 2, 12, 12, 5, 0, 7, 10],
                [14, 14, 18, 0, 10, 1, 10, 15, 7, 19],
                [9, 19, 5, 18, 7, 6, 14, 10, 15, 13],
                [18, 2, 1, 18, 10, 18, 4, 5, 13, 11],
                [1, 1, 2, 15, 2, 19, 17, 18, 16, 11],
                [15, 0, 8, 12, 7, 8, 9, 3, 5, 16]], aug=True, mod=mods)

    m.choose_pivots_mod()

    def test_for_error():
        b_count = 0
        while True:
            b = Matrix(rows=6, rand=True, mod=mods, aug=True)
            try:
                if b.choose_pivots_mod() is None:
                    if b_count == 0:
                        print("valid b to test")
                        print(b)
                    b_count += 1
            except:
                print("error happened")
                print(b)
                print(mods)
                break


if __name__ == '__main__':
    test_mod2()
