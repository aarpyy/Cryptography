from cryptography318.linear_algebra import *
from cryptography318.tools import *
from cryptography318.array import *
from random import randrange
import pytest
import timeit


def matmul(obj1, obj2, row_type=None, obj_type=None, mod=None):
    if row_type is None:
        row_type = list
    if obj_type is None:
        obj_type = list

    # attempts to transpose using built in methods, manually performs if no method exists
    transpose = getattr(obj2, 'transpose', None)
    if transpose is None:
        T = []
        for j in range(len(obj2[0])):
            T.append([])
            for i in range(len(obj2)):
                T[j].append(obj2[i][j])
    else:
        T = transpose()

    if isnumber(mod):
        if row_type == ArrayMod:
            return obj_type(map(lambda r: ArrayMod(map(lambda c: dot(c, r, mod), T), mod), obj1))
        return obj_type(map(lambda r: row_type(map(lambda c: dot(c, r, mod), T)), obj1))
    return obj_type(map(lambda r: row_type(map(lambda c: dot(c, r), T)), obj1))


def test_change_basis():
    S = LinearMap([[1, 0],
                   [0, -1]])
    B = Matrix([[1, -2],
                [2, 1]]).invert()
    T = Matrix([[-.6, .8],
                [.8, .6]])
    assert S.in_basis(B) == T


def test_ranknull(it=50):
    for _ in range(it):
        x = randrange(2, 10)
        m = Matrix(rows=x, cols=x, rand=True)
        assert m.rank() + m.null() == m.dimension()


def test_rref(it=50):
    for _ in range(it):
        x = randrange(2, 10)
        m = Matrix(rows=x, cols=x, rand=True)
        assert m.rref().is_rref()


def test_basis():
    A = LinearMap([[2, 2],
                   [1, 3]])
    assert A.is_eigen_value(1)


def test_linear_solve(it=50):
    for _ in range(it):
        m = Matrix(rand=True, aug=True)
        m = m.rref()
        if m.is_solvable():
            assert m.solve()


def test_inner_product():

    # works with two column vectors
    m = Matrix([[1, 2, 3]]).transpose()
    a = numpy.array([[1, 2, 3]]).transpose()
    r = m.inner(a)
    assert r == 14

    # works with one column vector and one row vector
    m = Matrix([[1, 2, 3]])
    a = numpy.array([[1, 2, 3]]).transpose()
    r = m.inner_prod(a)
    assert r == 14

    # works with two non-nested row vectors
    m = Matrix([1, 2, 3])
    a = numpy.array([1, 2, 3])
    r = m.inner_prod(a)
    assert r == 14

    matrix = Matrix([[1, 0, 0],
                     [0, 1, 0],
                     [0, 0, 1],
                     [-1, -1, -1]])
    t = matrix.T()
    assert Matrix.inner_prod(t[0]) == 2


def test_orthogonal():
    m = Matrix([[1, 1, 1, 1],
                [1, 1, -1, -1],
                [1, -1, 1, -1],
                [1, -1, -1, 1]])

    assert m.orthogonal()


def test_change():
    T = LinearMap([[4, -1, 1],
                   [3, -1, 1]])
    C = Matrix([[1, 0, 1],
                [3, 2, 0],
                [0, 3, -4]])


def test_remove_null_row():
    m = Matrix([[1, 2, 3],
                [4, 5, 6],
                [0, 0, 0],
                [0, 1, 0],
                [1, 2, 3],
                [0, 0, 0]])
    m.remove_null_row()
    assert len(m) == 4
    for row in m:
        assert any(aslist(row))


def test_remove_null_col():
    m = Matrix([[1, 0, 2, 3, 0, 0, 0],
                [5, 0, 4, 0, -0, 0, 2],
                [1, 0, 3, 0, 0, 0, 0]])
    m.remove_null_column()
    assert len(m[0]) == 4
    for column in m.T():
        assert any(aslist(column))


@pytest.mark.skip
def test_optimal_array_addition_mod():

    # testing which method is faster for adding arrays when there is a mod
    m = Matrix([[1, 2, 3, 4],
                [2, 3, 4, 5],
                [3, 4, 5, 6],
                [4, 5, 6, 7],
                [5, 6, 7, 8]])

    def sum_rows1(mat):
        array = numpy.array([0] * 4)
        for r in mat:
            array = (array + r) % 3

    def sum_rows2(mat):
        array = numpy.array([0] * 4)
        for r in mat:
            array += r
            if any(n > 3 for n in array):
                array %= 3

    def sum_rows3(mat):
        array = numpy.array([0] * 4)
        for r in mat:
            array += r
        mat.array %= 3

    mod_while = timeit.timeit(lambda: sum_rows1(m), number=100000)
    mod_if = timeit.timeit(lambda: sum_rows2(m), number=100000)
    mod_after = timeit.timeit(lambda: sum_rows3(m), number=100000)
    print(f"Modding during addition took {mod_while:.2f}s")
    print(f"Modding with if statement took {mod_if:.2f}s")
    print(f"Modding with after addition took {mod_after:.2f}s")


def test_norm():
    for _ in range(50):
        v = Matrix(rows=1, cols=3, rand=True)
        assert abs(pow(v.norm(), 2) - v.inner_prod(v)) < .5


def test_in_ortho_basis():
    T = LinearMap([[0, 1, 0],
                   [0, 0, 1],
                   [1, 0, 0]])
    B = Matrix([[1/sqrt(2), 1/sqrt(6)],
                [0, -2/sqrt(6)],
                [-1/sqrt(2), 1/sqrt(6)]])
    # print(T.in_ortho_basis(B))


def test_coordinates():
    M = Matrix([1, 2, -3]).transpose()
    B = Matrix([[1 / sqrt(2), 1 / sqrt(6)],
                [0, -2 / sqrt(6)],
                [-1 / sqrt(2), 1 / sqrt(6)]])
    r = M.coordinates(B)
    B = B.transpose()
    result = Matrix([0, 0, 0])
    for i in range(len(r)):
        result[0] += r[i][0] * B[i]

    assert result == M.T()


def test_to_fraction():
    M = Matrix([1, 2, -3]).transpose()
    B = Matrix([[1 / sqrt(2), 1 / sqrt(6)],
                [0, -2 / sqrt(6)],
                [-1 / sqrt(2), 1 / sqrt(6)]])
    r = M.coordinates(B)
    f = r.to_fraction()

    assert isinstance(f[0][0], str) and '/' in f[0][0]

    for i, row in enumerate(f):
        assert abs(evaluate(row[0]) - r[i][0]) < 0.1


def test_orthonormalize():
    matrix = Matrix([[1, 0, 0],
                     [0, 1, 0],
                     [0, 0, 1],
                     [-1, -1, -2]])
    k = matrix.orthonormalize(steps=True)


def test_eigen():
    A = Matrix([[-1, -3, 1],
                [3, 3, 1],
                [3, 0, 4]])
    A1 = A - Matrix(rows=3, cols=3, identity=1)
    A2 = A - Matrix(rows=3, cols=3, identity=2)
    A3 = A - Matrix(rows=3, cols=3, identity=3)

    v1 = Matrix([-1, 1, 1]).T()
    v2 = Matrix([-2, 3, 3]).T()
    v3 = Matrix([-3, 7, 9]).T()

    assert A * v1 == v1
    assert A * v2 == 2 * v2
    assert A * v3 == 3 * v3


def test_matmul():
    m = numpy.array([[2, 0, 5],
                     [1, 3, 4],
                     [9, -1, 4]])
    a = numpy.array([[1],
                     [3],
                     [-2]])
    print(matmul(m, a, list))
    print(Matrix(m) * Matrix(a))


if __name__ == '__main__':
    test_change_basis()
    test_ranknull()
    test_rref()
    test_basis()
    test_linear_solve()
    test_inner_product()
    test_orthogonal()
    test_change()
    test_remove_null_row()
    test_remove_null_col()
    # test_optimal_array_addition_mod()
    test_norm()
    test_in_ortho_basis()
    test_coordinates()
    test_to_fraction()
    # test_orthonormalize()
    test_eigen()
    test_matmul()
    m = Matrix([[1, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 1, 0, 0, 1]])
    c = ArrayMod([1, 2, 3], 5)
    a = Matrix([1, 2, 3]).T()
    b = [[1, 2, 3]]
    res = matmul(b, a, row_type=type(c), obj_type=Matrix, mod=5)
    print(res, type(res[0]))
