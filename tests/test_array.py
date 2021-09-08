from random import randrange
from unittest import TestCase

from cryptography318.array import *
from functools import reduce
import numpy
import pytest


def random_list(length=None):
    if length is None:
        length = randrange(3, 7)
    return [randrange(-15, 15) for _ in range(length)]


def list_equals(obj1, obj2):
    if isinstance(obj1[0], list):
        obj1 = obj1[0]
    if isinstance(obj2[0], list):
        obj2 = obj2[0]

    for i, e in enumerate(obj1):
        if e != obj2[i]:
            return False
    return True


def test_constructor():
    with pytest.raises(TypeError) as exc_info:
        Array('hello world')
    assert "requires input object" in str(exc_info.value)

    with pytest.raises(TypeError) as exc_info:
        Array([[1, 2, 3], [1, 2, 3]])
    assert "requires input object" in str(exc_info.value)

    assert Array(5).array == [5]
    assert Array(numpy.array([1, 2, 3])).array == [1, 2, 3]

    a = Array(numpy.int64(5))
    assert isinstance(a.array[0], int)
    a = Array([[1, 2, 3]])
    assert a == [1, 2, 3]
    a = Array([])
    assert a == []


def test_add():
    a = Array([1, 2, 3])
    b = Array([1, 2, 3])

    assert list_equals(a + b, Array([2, 4, 6]))

    for _ in range(50):
        length = randrange(3, 7)
        lst1 = random_list(length)
        lst2 = random_list(length)
        n1 = numpy.array(lst1)
        n2 = numpy.array(lst2)
        a1 = Array(lst1)
        a2 = Array(lst2)
        assert list_equals(n1 + n2, a1 + a2)


def test_str_repr():
    a = Array([1, 2, 3])
    assert str(a) == '[1, 2, 3]'
    assert repr(a) == 'Array([1, 2, 3])'


def test_eq():
    a = Array([1, 2, 3, 2, 4])
    b = Array([1, 2, 3, 2, 4])
    c = Array([1, 2, 3, 4, 2])

    assert a == b
    assert not a == c

    assert (a == 2) == [0, 1, 0, 1, 0]
    assert (a != 2) == [1, 0, 1, 0, 1]
    assert (a == {2, 3}) == [0, 1, 1, 1, 0]
    assert (a != {2, 3}) == [1, 0, 0, 0, 1]
    assert (2 in a)
    assert not (5 in a)

    for _ in range(50):
        length = randrange(3, 7)
        lst1 = random_list(length)
        lst2 = random_list(length)
        n1 = numpy.array(lst1)
        n2 = numpy.array(lst2)
        a1 = Array(lst1)
        a2 = Array(lst2)
        assert all(map(lambda e1, e2: bool(e1) == bool(e2), n1 < n2, a1 < a2))
        assert all(map(lambda e1, e2: bool(e1) == bool(e2), n1 <= n2, a1 <= a2))
        assert all(map(lambda e1, e2: bool(e1) == bool(e2), n1 > n2, a1 > a2))
        assert all(map(lambda e1, e2: bool(e1) == bool(e2), n1 >= n2, a1 >= a2))


def test_sub():
    a = Array([1, 2, 3])
    b = [1, 2, 4]

    assert list_equals(a - b, Array([0, 0, -1]))
    assert list_equals(b - a, Array([0, 0, 1]))

    for _ in range(50):
        length = randrange(3, 7)
        lst1 = random_list(length)
        lst2 = random_list(length)
        n1 = numpy.array(lst1)
        n2 = numpy.array(lst2)
        a1 = Array(lst1)
        a2 = Array(lst2)
        assert list_equals(n1 - n2, a1 - a2)


def test_mul():
    a = Array([1, 2, 3])
    b = [1, 2, 4]

    assert list_equals(a * b, Array([1, 4, 12]))
    assert list_equals(b * a, Array([1, 4, 12]))
    assert list_equals(a * 2, Array([2, 4, 6]))
    assert list_equals(2 * a, Array([2, 4, 6]))

    for _ in range(50):
        length = randrange(3, 7)
        lst1 = random_list(length)
        lst2 = random_list(length)
        n1 = numpy.array(lst1)
        n2 = numpy.array(lst2)
        a1 = Array(lst1)
        a2 = Array(lst2)
        assert list_equals(n1 * n2, a1 * a2)


def test_div():
    a = Array([2, 4, 6])
    b = [1, 2, 3]

    assert list_equals(a / b, Array([2, 2, 2]))
    assert list_equals(b / a, Array([.5, .5, .5]))
    assert list_equals(b // a, Array([0, 0, 0]))
    assert list_equals(a / 2, Array([1, 2, 3]))
    assert list_equals(a / 3, Array([2 / 3, 4 / 3, 2]))
    assert list_equals(a // 3, Array([0, 1, 2]))

    for _ in range(50):
        length = randrange(3, 7)
        lst1 = random_list(length)
        lst2 = random_list(length)
        if 0 in lst2:
            continue
        n1 = numpy.array(lst1)
        n2 = numpy.array(lst2)
        a1 = Array(lst1)
        a2 = Array(lst2)
        assert list_equals(n1 / n2, a1 / a2)
        assert list_equals(n1 // n2, a1 // a2)

    a = ArrayMod([2, 4, 6, 8, 10], 12)
    assert (k := a // 2) == [1, 2, 3, 4, 5] and k.mod == 12
    assert (k := a / 2) == [1, 2, 3, 4, 5] and k.mod == 6

    with pytest.raises(ValueError) as exc_info:
        a // 4
    assert "does not divide" in str(exc_info.value)

    a = ArrayMod([2, 4, 6, 8, 10], 14)

    # can divide because row divisible by 2, then whats left to divide by is 3 which is invertible mod 14, so
    # whole process looks like dividing row by 2, then multiplying row by pow(3, -1, 14)
    assert (k := a // 6) == [5, 10, 1, 6, 11] and k.mod == 14

    a = ArrayMod([2, 4, 6, 8, 10], 11)
    assert (k := a // 2) == [1, 2, 3, 4, 5] and k.mod == 11

    with pytest.raises(ValueError) as exc_info:
        a / 2
    assert "not divisible" in str(exc_info.value)

    assert (k := a // 3) == [8, 5, 2, 10, 7] and k.mod == 11


def test_pow_mod():
    a = Array([1, 2, 3])

    assert list_equals(pow(a, 2), [1, 4, 9])
    assert list_equals(pow(a, 2, 5), [1, 4, 4])
    assert list_equals(pow(a, 6, 7), [1, 1, 1])

    a *= 5
    assert list_equals(a % 5, [0, 0, 0])
    assert list_equals(a % 4, [1, 2, 3])


def test_make_pivot():
    a = Array([1, 2, 3, 4])
    b = Array([4, 3, 2, 1])

    assert a.make_pivot(copy=True) == [1, 2, 3, 4]
    assert a.make_pivot(index=1, copy=True) == [.5, 1, 1.5, 2]
    assert b.make_pivot(copy=True) == [1, .75, .5, .25]

    c = Array([4, 1, 5, 1, 2, 9])
    # print(c.make_pivot(2, copy=True))


def test_where():
    a = Array([1, 2, 3, 1])

    assert where(a)[0] == [0, 1, 2, 3]
    assert where(a == 1)[0] == [0, 3]


def test_copy():
    a = Array([1, 2, 3])
    b = a.copy()
    assert a is not b and a == b

    c = a[:]
    assert isinstance(c, Array) and c is not a and c == a
    c[0] = 5
    assert c[0] == 5 and a[0] != 5

    a = Array([1, 2, 3])
    b = Array([4, 5, 6])
    a, b = b[:], a[:]
    assert a == [4, 5, 6] and b == [1, 2, 3]
    assert isinstance(a, Array) and isinstance(b, Array)

    a = Array([1, 2, 3])
    b = Array([4, 5, 6])
    a, b = b.copy(), a.copy()
    assert a == [4, 5, 6] and b == [1, 2, 3]


def test_shift():
    a = Array([1, 2, 3, 4])
    b = a.shift_elements(shift=2, copy=True)
    assert b == [3, 4, 1, 2]
    b.shift_elements(shift=1)
    assert a == [1, 2, 3, 4] and b == [2, 3, 4, 1]


if __name__ == '__main__':
    test_add()
    test_constructor()
    test_str_repr()
    test_eq()
    test_sub()
    test_mul()
    test_div()
    test_pow_mod()
    test_make_pivot()
    test_where()
    test_copy()
    test_shift()
    test_shift()
    a = ArrayMod([2, 4, 6, 8, 12], mod=14)
    g = gcd(*a.array)
    pivot = lambda arr: reduce(lambda i, c: c if (
                    gcd(c, arr.mod) == 1 or (g > 1 and gcd(c//g, arr.mod) == 1)
            ) else i, arr.array[::-1], False)
    print(p := pivot(a))
    print(a // p)
