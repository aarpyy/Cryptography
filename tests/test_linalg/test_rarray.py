from cryptography318.linalg.arrayabc import *
from typing import Iterator
from fractions import Fraction
from timeit import timeit
from random import randrange


def test_constructor():
    a = rarray([1, 2, 3, 4])
    assert isinstance(a, rarray) and isinstance(a, ABCArray)
    assert hasattr(a, '_array') and hasattr(a, '_dtype')


def test_dtype():
    a = rarray([1, 2, 3, 4])
    assert a.dtype is int
    b = a.astype(float)
    assert b.dtype is float and all(isinstance(e, float) for e in b.array)
    assert a.dtype is int and all(isinstance(e, int) for e in a.array)
    a.astype(float, _update=True)
    assert a.dtype is float and all(isinstance(e, float) for e in a.array)


def test_str():
    a = rarray([1, 2, 3, 4])
    assert str(a) == str(a.array)
    assert isinstance(eval(repr(a)), rarray)


def test_iter():
    a = rarray([1, 2, 3, 4])
    _iter = iter(a)
    assert isinstance(_iter, Iterator)
    assert _iter.__next__() == 1
    assert _iter.__next__() == 2
    assert _iter.__next__() == 3
    assert _iter.__next__() == 4


def test_set_item():
    a = rarray([1, 2, 3, 4])
    a[0] = 5
    assert a[0] == 5 and a[1:] == [2, 3, 4]

    a[0] = 5.0
    # float higher in TRO than int
    assert isinstance(a[0], float) and a == [5, 2, 3, 4] and a.dtype is float

    a[0] = complex(1, 1)
    # complex higher than float
    assert a[0] == complex(1, 1) and a.dtype is complex
    assert all(isinstance(e, complex) for e in a.array)

    a[0] = Fraction(1, 2)
    # complex > Fraction, so convert Fraction to complex
    assert a[0] == complex(.5, 0) and a.dtype is complex
    assert all(isinstance(e, complex) for e in a.array)

    a = rarray([1, 2, 3, 4])

    # Fraction > int
    a[0] = Fraction(1, 2)
    assert a[0] == Fraction(1, 2) and a.dtype is Fraction
    assert all(isinstance(e, Fraction) for e in a.array)

    a[0] = 1.0
    # float > Fraction
    assert a[0] == 1 and a.dtype is float
    assert all(isinstance(e, float) for e in a.array)


def test_add():
    a = rarray([1, 2, 3, 4])
    assert a + 2 == [3, 4, 5, 6]

    assert a + [1, 2, 3, 4] == [2, 4, 6, 8]

    b = rarray([-1, -2, -3, -4])
    assert a + b == [0, 0, 0, 0]

    c = [1, 2, 3, 4] + a
    assert c == [2, 4, 6, 8] and isinstance(c, rarray)

    c = a - 2.0
    assert c.dtype is float and c == [-1, 0, 1, 2]

    c = [1, 2, 3, 4] - a
    assert c == [0, 0, 0, 0] and isinstance(c, rarray)

    assert a - [1, 1, 1, 1] == [0, 1, 2, 3]


def test_mul():
    a = rarray([1, 2, 3, 4])
    assert a * 2 == [2, 4, 6, 8]
    assert 2 * a == [2, 4, 6, 8]

    assert a * [1, 2, 3, 4] == [1, 4, 9, 16]
    assert [0, 0, 0, 0] * a == [0, 0, 0, 0]


def test_div():
    a = rarray([2, 4, 6, 8])
    b = a // 2
    assert b == [1, 2, 3, 4] and isinstance(b, rarray) and b.dtype is int

    b = a / 2
    assert b == [1, 2, 3, 4] and isinstance(b, rarray) and b.dtype is float

    b = a // 3
    assert b == [0, 1, 2, 2] and isinstance(b, rarray) and b.dtype is int

    b = a / 3
    assert b == [2/3, 4/3, 2, 8/3] and isinstance(b, rarray) and b.dtype is float

    b = a // 2.0
    assert b == [1, 2, 3, 4] and b.dtype is float


def test_all():
    test_constructor()
    test_dtype()
    test_str()
    test_iter()
    test_set_item()
    test_add()
    test_mul()


if __name__ == '__main__':
    test_div()
