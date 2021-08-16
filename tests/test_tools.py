from cryptography318.tools import *


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


if __name__ == '__main__':
    test_replace_all()
    test_string_reduce()
