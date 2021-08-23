from random import randrange
from cryptography318.array import *
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
        ArrayMod('hello world', 5)
    assert "requires input object" in str(exc_info.value)

    with pytest.raises(TypeError) as exc_info:
        ArrayMod([], 5)
    assert "requires input object" in str(exc_info.value)

    with pytest.raises(TypeError) as exc_info:
        ArrayMod([[1, 2, 3]], 5)
    assert "requires input object" in str(exc_info.value)

    with pytest.raises(TypeError) as exc_info:
        ArrayMod([1, 2, 3])
    assert "required positional argument: \'mod\'" in str(exc_info.value)

    assert list_equals(ArrayMod([1, 2, 3], 4), [1, 2, 3])
    assert list_equals(ArrayMod([1, 2, 5], 3), [1, 2, 2])

    a = ArrayMod(numpy.int64(5), 4)
    assert a[0] == 1


def test_add():
    a = ArrayMod([1, 2, 3], 5)
    b = ArrayMod([1, 2, 3], 5)

    assert list_equals(a + b, [2, 4, 1])

    for _ in range(50):
        length = randrange(3, 7)
        lst1 = random_list(length)
        lst2 = random_list(length)
        n1 = numpy.array(lst1)
        n2 = numpy.array(lst2)
        a1 = ArrayMod(lst1, 5)
        a2 = ArrayMod(lst2, 5)
        assert list_equals((n1 + n2) % 5, a1 + a2)


def test_str_repr():
    a = ArrayMod([1, 2, 3], 5)
    assert str(a) == '[1, 2, 3]'
    assert repr(a) == 'ArrayMod([1, 2, 3], mod=5)'


def test_eq():
    a = ArrayMod([1, 2, 3, 2, 4], 5)
    b = ArrayMod([1, 2, 3, 2, 4], 5)
    c = ArrayMod([1, 2, 3, 4, 2], 5)
    d = ArrayMod([1, 2, 3, 4, 2], 6)

    assert a == b
    assert not a == c

    # finding invertible elements
    assert (d == 'inv') == [1, 0, 0, 0, 0]
    assert (a == 'inv') == [1, 1, 1, 1, 1]

    # arrays with different moduli
    assert not c == d

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
        n1 = numpy.array(lst1) % 5
        n2 = numpy.array(lst2) % 5
        a1 = ArrayMod(lst1, 5)
        a2 = ArrayMod(lst2, 5)
        assert all(map(lambda e1, e2: bool(e1) == bool(e2), n1 < n2, a1 < a2))
        assert all(map(lambda e1, e2: bool(e1) == bool(e2), n1 <= n2, a1 <= a2))
        assert all(map(lambda e1, e2: bool(e1) == bool(e2), n1 > n2, a1 > a2))
        assert all(map(lambda e1, e2: bool(e1) == bool(e2), n1 >= n2, a1 >= a2))


def test_sub():
    a = ArrayMod([1, 2, 3], 5)
    b = [1, 2, 4]

    assert list_equals(a - b, [0, 0, 4])
    assert list_equals(b - a, [0, 0, 1])

    for _ in range(50):
        length = randrange(3, 7)
        lst1 = random_list(length)
        lst2 = random_list(length)
        n1 = numpy.array(lst1)
        n2 = numpy.array(lst2)
        a1 = ArrayMod(lst1, 5)
        a2 = ArrayMod(lst2, 5)
        assert list_equals((n1 - n2) % 5, a1 - a2)


def test_mul():
    a = ArrayMod([1, 2, 3], 5)
    b = [1, 2, 4]

    assert list_equals(a * b, [1, 4, 2])
    assert list_equals(b * a, [1, 4, 2])
    assert list_equals(a * 2, [2, 4, 1])
    assert list_equals(2 * a, [2, 4, 1])

    for _ in range(50):
        length = randrange(3, 7)
        lst1 = random_list(length)
        lst2 = random_list(length)
        n1 = numpy.array(lst1)
        n2 = numpy.array(lst2)
        a1 = ArrayMod(lst1, 5)
        a2 = ArrayMod(lst2, 5)
        assert list_equals((n1 * n2) % 5, a1 * a2)


def test_div():
    a = ArrayMod([2, 4, 6], 5)

    assert list_equals(a / 2, [1, 2, 3])
    assert list_equals(a // 2, [1, 2, 3])

    b = ArrayMod([3, 6, 6], 9)
    k = b / 6
    assert list_equals(k, [2, 1, 1]) and k.mod == 3


def test_pow_mod():
    a = ArrayMod([1, 2, 3], 11)

    assert list_equals(pow(a, 2), [1, 4, 9])
    assert list_equals(pow(a, 2, 5), [1, 4, 4])
    assert list_equals(pow(a, -1), [1, 6, 4])

    a *= 5
    assert list_equals(a, [5, 10, 4])


def test_make_pivot():
    a = ArrayMod([1, 2, 3, 4], 7)

    assert a.make_pivot(index=2)[2] == 1

    assert a.make_pivot()[0] == 1


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
