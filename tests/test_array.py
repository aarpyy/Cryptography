import pytest
from cryptography318.linear_algebra import *


def test_constructor():
    # create object
    m = array_mod([2, 3, 5], mod=5)
    assert isinstance(m, array_mod)

    # error thrown with no given mod
    with pytest.raises(ValueError) as exc_info:
        array_mod([1, 2, 3])
    assert "modulus" in str(exc_info.value)

    # error thrown with multi-row input
    with pytest.raises(ValueError) as exc_info:
        array_mod([[1, 2], [3, 4]], mod=5)
    assert "row vector" in str(exc_info.value)

    # inheritance of modulus
    a = array_mod(m)
    assert a.mod == m.mod

    # construction without given array
    a = array_mod(cols=3, mod=4)
    assert isinstance(a, array_mod)

    # construction from numpy.ndarray
    a = array_mod(numpy.array([[1, 2, 3]]), mod=4)
    assert isinstance(a, array_mod)

    # construction from Matrix
    a = array_mod(Matrix([[1, 2, 3]]), mod=4)
    assert isinstance(a, array_mod)


def test_get_item():
    m = array_mod([2, 3, 5], mod=5)
    assert m[0] == 2
    with pytest.raises(IndexError):
        var = m[3]


def test_set_item():
    m = array_mod([2, 3, 5], mod=5)
    assert m[0] == 2
    m[0] = 4
    assert m[0] == 4

    # can only set number values
    with pytest.raises(ValueError) as exc_info:
        m[0] = 'h'
    assert "<class 'str'> unacceptable" in str(exc_info.value)

    # float values get converted to integers
    m[0] = 3.1
    assert m[0] == 3


def test_len():
    m = array_mod([2, 3, 5], mod=5)
    assert len(m) == 3


def test_iter():
    m = array_mod([2, 3, 5], mod=5)
    for n in m:
        assert isinstance(n, int)

    # object returned by list is actual list, not layered array_mod inside list
    assert isinstance(list(m), list)


def test_str():
    m = array_mod([2, 3, 5], mod=5)
    assert isinstance(str(m), str)

    assert '2' in str(m)


def test_add():
    m = array_mod([2, 3, 5], mod=5)
    a = array_mod([7, 5, 2], mod=6)
    b = a + m
    c = m + a

    # test if values are as expected
    for e in b:
        assert e in [1, 2, 3]
    for e in c:
        assert e in [2, 3, 4]

    # test if modulus was inherited correctly
    assert b.mod == a.mod
    assert c.mod == m.mod

    # adding number to array
    b = m + 4
    for e in b:
        assert e in [1, 2, 4]

    # prevent adding array to number
    with pytest.raises(ValueError) as exc_info:
        b = 4 + m
    assert "Unable to add" in str(exc_info.value)

    # prevent adding multi-dimensional array
    with pytest.raises(AttributeError) as exc_info:
        m += numpy.array([[1, 2], [3, 4]])
    assert "multi-dimensional array is unsupported" in str(exc_info.value)

