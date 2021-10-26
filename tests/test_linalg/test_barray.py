from cryptography318.linalg.arrayabc import *
from timeit import timeit
from random import randrange
from math import prod


def test_invert():
    a = bitarray(randrange(32), _len=5)
    b = a.complement()
    assert a + b == [1, 1, 1, 1, 1]
    assert all(x != y for x, y in zip(a, b))


def test_as_integer():
    assert bitarray.as_integer([1, 0, 0, 1]) == 9
    assert bitarray.as_integer([1, 1, 0, 1]) == 13


def test_repr():
    a = bitarray(randrange(32), _len=5)
    b = eval(repr(a))
    assert b == a and isinstance(b, bitarray)


def test_iter():
    a = bitarray(randrange(32), _len=5)
    for bit in a:
        assert bit in (0, 1)


def test_set_item():
    a = bitarray(6, _len=8)
    assert a == [0, 0, 0, 0, 0, 1, 1, 0]
    a[1] = 1
    assert a == [0, 1, 0, 0, 0, 1, 1, 0]
    a[0] = 1
    assert a == [1, 1, 0, 0, 0, 1, 1, 0]
    a[6] = 0
    assert a == [1, 1, 0, 0, 0, 1, 0, 0]
    a[6] = 0
    assert a == [1, 1, 0, 0, 0, 1, 0, 0]
    a[1] = 3
    assert a == [1, 1, 0, 0, 0, 1, 0, 0]


# __add__ == __radd__ == __sub__ == __rsub__
def test_add():
    a = bitarray(randrange(32), _len=5)
    b = bitarray(randrange(32), _len=5)
    c = a + b
    for i, bit in enumerate(c):
        assert (a[i] + b[i]) & 1 == bit


def test_mul():
    a = bitarray(13, _len=5)
    b = bitarray(21, _len=5)
    # dot product of 01101 and 10101 = 1*1 + 1*1
    assert a * b == [0, 0, 1, 0, 1]

    a = bitarray(17, _len=5)
    b = bitarray(30, _len=5)
    assert a * b == [1, 0, 0, 0, 0]


def test_pos_neg():
    a = bitarray(17, _len=5)
    assert +a == a == -a


def test_append():
    a = bitarray(randrange(32), _len=5)
    b = a[:]
    assert len(a) == 5
    a.append(1)
    assert len(a) == 6
    assert a == list(b) + [1]


def test_reverse():
    a = bitarray(29, _len=6)
    a.reverse()
    assert a == [1, 0, 1, 1, 1, 0]
    a.reverse()
    assert a == [0, 1, 1, 1, 0, 1]


def test_extend():
    a = bitarray(17, _len=5)
    a.extend([1, 1, 0, 1])
    assert a == [1, 0, 0, 0, 1, 1, 1, 0, 1]
    assert len(a) == 9


def test_pop():
    a = bitarray(21, _len=5)
    assert a.pop()
    assert a == [0, 1, 0, 1] and len(a) == 4
    assert not a.pop(2)
    assert a == [0, 1, 1] and len(a) == 3


def test_remove():
    a = bitarray(21, _len=5)
    a.remove(1)
    assert a == [0, 1, 0, 1] and len(a) == 4

    a.remove(1)
    assert a == [0, 0, 1] and len(a) == 3

    a.remove(0)
    assert a == [0, 1] and len(a) == 2


def test_insert():
    a = bitarray(17, _len=5)
    a.insert(2, 1)
    assert a == [1, 0, 1, 0, 0, 1] and len(a) == 6

    a.insert(0, 1)
    assert a == [1, 1, 0, 1, 0, 0, 1] and len(a) == 7

    a.insert(7, 0)
    assert a == [1, 1, 0, 1, 0, 0, 1, 0] and len(a) == 8


def test_get_item():
    a = bitarray(21, _len=5)
    assert 1 == a[0] == a[2] == a[4]
    assert 0 == a[1] == a[3]


def test_compare_list():

    def dot(a, b):
        return reduce(lambda r, c: r + prod(c), zip(a, b), 0)

    time_b, time_l = 0, 0
    for _ in range(pow(10, 3)):
        a = bitarray(randrange(32), _len=25)
        b = bitarray(randrange(32), _len=25)
        x = a.array
        y = b.array
        time_b += timeit(lambda: a.dot(b), number=pow(10, 3))
        time_l += timeit(lambda: dot(x, y), number=pow(10, 3))
    print(f"barray: {time_b:.2f}s")
    print(f"list: {time_l:.2f}s")


def test_all():
    test_append()
    test_add()
    test_repr()
    test_iter()
    test_invert()
    test_set_item()
    test_as_integer()
    test_mul()
    test_pos_neg()
    test_reverse()
    test_extend()
    test_pop()
    test_remove()
    test_insert()
    test_get_item()


if __name__ == '__main__':
    test_all()
