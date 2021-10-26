from cryptography318.linalg.array_old import *
from cryptography318.linalg.linear_algebra import Matrix as Matrix2
import operator
from cryptography318.core.tools import *
import numpy
import timeit
from random import choice


def list_equals(obj1, obj2):
    for i, e in enumerate(obj1):
        if e != obj2[i]:
            return False
    return True


def test_add():
    b = Matrix([[1, 2, 3],
                [3, 2, 1]])

    assert isinstance(b, Matrix) and isinstance(b[0], Array)

    c = Matrix([[1, 2, 3],
                [1, 0, 1]])

    a = b + c
    assert isinstance(a, Matrix) and isinstance(a[0], Array)
    assert a == [[2, 4, 6], [4, 2, 2]]

    a += 5
    assert a == [[7, 9, 11], [9, 7, 7]]


def test_sub():
    b = Matrix([[1, 2, 3],
                [3, 2, 1]])
    c = Matrix([[1, 2, 3],
                [1, 0, 1]])
    a = b - c
    assert isinstance(a, Matrix) and isinstance(a[0], Array)

    assert a == [[0, 0, 0], [2, 2, 0]]

    a -= 2
    assert a == [[-2, -2, -2], [0, 0, -2]]


def test_random_op():
    ops_sq = [operator.add, operator.sub]
    for _ in range(50):
        cols = randrange(2, 10)
        a = Matrix(cols=cols, rand=True)
        b = Matrix(rows=cols, rand=True)
        a2 = Matrix2(a.tolist())
        b2 = Matrix2(b.tolist())
        scalar = randrange(2, 50)
        if len(a) == len(a[0]) and len(b) == len(b[0]):
            op = choice(ops_sq)
            c = op(a, b) + scalar
            c2 = op(a2, b2) + scalar
        else:
            c = a * b
            c2 = a2 * b2
        dt = a / scalar
        df = a // scalar
        dt2 = a2 / scalar
        df2 = a2 // scalar
        assert len(c) == len(a) and len(c[0]) == len(b[0])
        assert len(c2) == len(c) and len(c2[0]) == len(c[0])
        assert c == c2
        assert dt == dt2 and df == df2


def test_inverse(time=False):
    for _ in range(50):
        length = randrange(2, 8)
        a = Matrix(rows=length, cols=length, rand=True)
        b = Matrix2(a.tolist())
        assert a.complement() == b.invert()
        assert isinstance(a, Matrix) and isinstance(a[0], Array) and isinstance(a[0][0], (int, float))

    if time:
        mat = timeit.timeit(lambda: a.complement(), number=10000)
        mat2 = timeit.timeit(lambda: b.invert(), number=10000)
        print(f"Matrix took: {mat:.2f}s")
        print(f"Matrix2 took: {mat2:.2f}s")


def test_det(time=False):
    for _ in range(50):
        length = randrange(2, 8)
        a = Matrix(rows=length, cols=length, rand=True)
        b = numpy.array(a.tolist())
        assert abs(a.det() - numpy.linalg.det(b)) < .5

    if time:
        length = randrange(2, 8)
        a = Matrix(rows=length, cols=length, rand=True)
        b = numpy.array(a.tolist())
        mat = timeit.timeit(lambda: a.det(), number=10000)
        np = timeit.timeit(lambda: numpy.linalg.det(b), number=10000)
        print(f"matrix took: {mat:.2f}s")
        print(f"numpy took: {np:.2f}s")


def test_append():
    a = Matrix([[1, 2, 3],
                [4, 5, 6]])
    a.append([4, 6, 7])
    assert a == [[1, 2, 3], [4, 5, 6], [4, 6, 7]]


def test_copy():
    a = Matrix([[1, 2, 3],
                [4, 5, 6]])
    b = a.copy()


def test_cross():
    for _ in range(50):
        length = randrange(2, 8)
        a = Matrix(rows=length - 1, cols=length, rand=True)
        b = a.cross()
        for e in a:
            assert not dot(e, b)


def test_T():
    for _ in range(50):
        rows = randrange(3, 8)
        cols = randrange(2, 8)
        a = Matrix(rows=rows, cols=cols, rand=True)
        assert list_equals(transpose_obj(a), a.T)
        a.remove_row(1)
        assert list_equals(transpose_obj(a), a.T)
        a.remove_column(1)
        assert list_equals(transpose_obj(a), a.T)

        a = Matrix(rows=rows, cols=cols, rand=True)
        assert list_equals(transpose_obj(a), a.T)
        a = a.remove_row(1, copy=True)
        assert list_equals(transpose_obj(a), a.T)
        a = a.remove_column(1, copy=True)
        assert list_equals(transpose_obj(a), a.T)

        b = Matrix(rows=rows, cols=rows, identity=1)
        b[1][1] = 0
        b[1][2] = 1
        b[2][2] = 0
        b.remove_null_column()
        assert list_equals(transpose_obj(b), b.T)

        b.remove_null_row()
        assert list_equals(transpose_obj(b), b.T)

        b = Matrix(rows=rows, cols=rows, identity=1)
        b[1][1] = 0
        b[1][2] = 1
        b[2][2] = 0
        b = b.remove_null_column(copy=True)
        assert list_equals(transpose_obj(b), b.T)

        b = b.remove_null_row(copy=True)
        assert list_equals(transpose_obj(b), b.T)


if __name__ == '__main__':
    # test_add()
    # test_random_op()
    # test_det()
    # test_append()
    # test_copy()
    # test_cross()
    # test_inverse()
    x = Symbol('x')
    b = Matrix([[0, 2, 1],
                [-1, -2, 0],
                [1, 1, -1]])
    print((b - Matrix(rows=3, identity=x)).minor())
    print(b.eigvals())
