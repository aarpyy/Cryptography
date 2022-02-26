from cryptography318.core.tools import *
from cryptography318.deprecated.linear_algebra import *


def test_replace_all():
    s = 'abcde;jkl\'\nghk'
    to_replace = ';\'\n'
    assert all(c not in to_replace for c in replace_all(s, to_replace))


def test_string_reduce():
    n = 2.0001
    assert '2' == string_reduce(n)

    n = 2.000
    assert '2' == string_reduce(n)

    n = -2.001
    assert '-2.001' == string_reduce(n)


def test_evaluate():
    M = Matrix([1, 2, -3]).transpose()
    B = Matrix([[1 / sqrt(2), 1 / sqrt(6)],
                [0, -2 / sqrt(6)],
                [-1 / sqrt(2), 1 / sqrt(6)]])
    r = M.coordinates(B)
    f = r.to_fraction()

    # evaluate can eval both matrices and string-fractions
    assert evaluate(f) == r
    assert evaluate(f[0][0]) == r[0][0]
    assert evaluate(f[1][0]) == r[1][0]


if __name__ == '__main__':
    test_replace_all()
    test_string_reduce()
    test_evaluate()
