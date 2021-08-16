from cryptography318.linear_algebra import *
from random import randrange
import pytest
import timeit


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
    r = m.inner_prod(a)
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
