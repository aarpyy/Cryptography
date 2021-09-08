from cryptography318.fraction import *
import pytest
import fractions
from random import randrange
from timeit import timeit
from sympy import factorint as syfactor
from cryptography318.linalg import matrix


def test_constructor():
    x = SquareRoot(2)
    assert isinstance(x, SquareRoot) and x.value == sqrt(2)

    # construct from complex that is int
    x = SquareRoot(complex(2, 0))
    assert isinstance(x, SquareRoot) and x.value == sqrt(2)

    # construct from float that is int
    x = SquareRoot(float(2.0))
    assert isinstance(x, SquareRoot) and x.value == sqrt(2)

    # construct from fractions.Fraction object that is int
    x = SquareRoot(fractions.Fraction(2, 1))
    assert isinstance(x, SquareRoot) and x.value == sqrt(2)

    # non-integer inside radical invalid, float of approximate value should be returned
    x = SquareRoot(2.1, 3.0)
    assert isinstance(x, float) and x == 3 * sqrt(2.1)

    x = SquareRoot('sqrt(12)')
    assert x.value == sqrt(12) and x._value == 3

    x = SquareRoot('sqrt(8)')
    assert x._value == 2 and x._scalar == 2

    # cannot create sqrt object from incomplete string
    with pytest.raises(ValueError) as exc_info:
        x = SquareRoot('sqrt(8')
    assert "invalid string literal" in str(exc_info.value)

    # reduces to integer value of 4
    x = SquareRoot('2*sqrt(4)')
    assert x == 4 and isinstance(x, int)

    x = SquareRoot('2 * sqrt(4.0)')
    assert x == 4 and isinstance(x, int)

    # matches even with whitespace, converts floats to ints
    x = SquareRoot('   2.0    *   sqrt(4.0)  ')
    assert x == 4 and isinstance(x, int)

    x = SquareRoot('2.1 * sqrt(4.0)')
    assert x == 4.2

    # can't have decimal value with leading 0. as radical value
    with pytest.raises(ValueError) as exc_info:
        x = SquareRoot('.1 * sqrt(.2)')
    assert "invalid string literal" in str(exc_info.value)

    # can have decimal value with leading 0. as coefficient
    x = SquareRoot('.1 * sqrt(2)')
    assert x == sqrt(2) * .1


def test_add():
    x = SquareRoot(2)
    assert x + 2 == sqrt(2) + 2

    assert x + SquareRoot(2, 5) == SquareRoot(2, 6)


def test_mul():
    x = SquareRoot(2)
    assert x * SquareRoot(2) == 2
    assert x * SquareRoot(4) == SquareRoot(8)
    assert str(x * 2) == '(2*sqrt(2))'

    assert x * SquareRoot(4, 2) == 4 * SquareRoot(2)
    assert SquareRoot(2, 2) * SquareRoot(4, 2) == 8 * SquareRoot(2)

    assert SquareRoot(5, 2) * SquareRoot(10, 2) == SquareRoot(2, 20)


def test_div():
    x = SquareRoot(22, 2)
    assert x / SquareRoot(11, 4) == 0.5 * SquareRoot(2)

    assert x // SquareRoot(33) == 1
    assert x / SquareRoot(33) == 2 * sqrt(2/3)


def reduce_sqrt1(operand):
    s = 1
    factors = factor(operand)

    if factors is None or factors == {operand: 1}:
        return operand, s

    for f in factors:
        if factors[f] > 1:
            s *= pow(f, factors[f] // 2)
            operand //= pow(f, factors[f])
    return operand, s


def reduce_sqrt2(operand):
    s = 1
    factors = syfactor(operand)

    if factors is None or factors == {operand: 1}:
        return operand, s

    for f in factors:
        if factors[f] > 1:  # if there is a square or greater in the factors, push it to the coeff
            s *= pow(f, factors[f] // 2)
            operand //= pow(f, factors[f])
    return operand, s


def test_reduce_sq():
    time_1, time_2 = 0, 0
    for _ in range(1000):
        n = randrange(5, 5000)
        time_1 += timeit(lambda: reduce_sqrt1(n), number=1000)
        time_2 += timeit(lambda: reduce_sqrt2(n), number=1000)
    print(f"time for-loop1: {time_1:.2f}s")
    print(f"time for-loop2: {time_2:.2f}s")


def test_constructor_timeit():
    time_sr, time_fr = 0, 0
    for _ in range(10000):
        v1, v2 = randrange(5, 50), randrange(5, 50)
        time_sr += timeit(lambda: SquareRoot(v1, v2), number=1000)
        time_fr += timeit(lambda: fractions.Fraction(v1, v2), number=1000)

    print(f"time sr: {time_sr:.2f}s")
    print(f"time fr: {time_fr:.2f}s")


if __name__ == '__main__':
    test_constructor()
    test_div()
    test_mul()
    test_add()
    test_reduce_sq()

