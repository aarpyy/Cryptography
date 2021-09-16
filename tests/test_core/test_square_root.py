from cryptography318.core.expr import Sqrt
import pytest
import fractions
from math import sqrt
from random import randrange
from timeit import timeit
from cryptography318.numbers.factor import factor


def test_constructor():
    x = Sqrt(2)
    assert isinstance(x, Sqrt) and x.eval == sqrt(2)

    # construct from float that is int
    x = Sqrt(float(2.0))
    assert isinstance(x, Sqrt) and x.eval == sqrt(2)

    # construct from fractions.Fraction object that is int
    x = Sqrt(fractions.Fraction(2, 1))
    assert isinstance(x, Sqrt) and x.eval == sqrt(2)

    # non-integer inside radical invalid, float of approximate value should be returned
    x = Sqrt(2.1, 3.0)
    assert isinstance(x, float) and x == 3 * sqrt(2.1)

    x = Sqrt('sqrt(12)')
    assert x.eval == sqrt(12) and x.radicand == 3

    x = Sqrt('sqrt(8)')
    assert x.radicand == 2 and x.integer == 2

    # cannot create sqrt object from incomplete string
    with pytest.raises(ValueError) as exc_info:
        x = Sqrt('sqrt(8')
    assert "invalid string literal" in str(exc_info.value)

    # reduces to integer value of 4
    x = Sqrt('2*sqrt(4)')
    assert x == 4 and isinstance(x, Sqrt)
    assert x.radicand == 1 and x.integer == 4

    x = Sqrt('2 * sqrt(4.0)')
    assert x == 4 and isinstance(x, Sqrt)
    assert x.radicand == 1 and x.integer == 4

    # matches even with whitespace, converts floats to ints
    x = Sqrt('   2.0    *   sqrt(4.0)  ')
    assert x == 4 and isinstance(x, Sqrt)
    assert x.radicand == 1 and x.integer == 4

    x = Sqrt('2.1 * sqrt(4.0)')
    assert x == 4.2

    # can't have decimal value with leading 0. as radical value
    with pytest.raises(ValueError) as exc_info:
        x = Sqrt('.1 * sqrt(.2)')
    assert "invalid string literal" in str(exc_info.value)

    # can have decimal value with leading 0. as coefficient
    x = Sqrt('.1 * sqrt(2)')
    assert x == sqrt(2) * .1


def test_add():
    x = Sqrt(2)
    assert x + 2 == sqrt(2) + 2

    assert x + Sqrt(2, 5) == Sqrt(2, 6)


def test_mul():
    x = Sqrt(2)
    assert x * Sqrt(2) == 2
    assert x * Sqrt(4) == Sqrt(8)
    assert str(x * 2) == '(2*sqrt(2))'

    assert x * Sqrt(4, 2) == 4 * Sqrt(2)
    assert Sqrt(2, 2) * Sqrt(4, 2) == 8 * Sqrt(2)

    assert Sqrt(5, 2) * Sqrt(10, 2) == Sqrt(2, 20)


def test_div():
    x = Sqrt(22, 2)
    assert x / Sqrt(11, 4) == 0.5 * Sqrt(2)

    assert x // Sqrt(33) == 1
    assert x / Sqrt(33) == 2 * sqrt(2 / 3)

    assert str(10 // Sqrt(5, 2)) == 'sqrt(5)'

    assert str(Sqrt(5, 20) // Sqrt(10, 2)) == '(5*sqrt(2))'

    a = Sqrt(28)
    b = 56 / a
    assert isinstance(b, Sqrt) and b == Sqrt(7, 4)


def test_all():
    test_div()
    test_mul()
    test_add()
    test_constructor()


def reduce_sqrt1(operand, s=1):
    factors = factor(operand)

    if factors == {operand: 1} or not factors:
        return operand, s

    for f in factors:
        if factors[f] > 1:  # if there is a square or greater in the factors, push it to the coeff
            power = pow(f, factors[f] // 2)
            s *= power
            operand //= (power * power)
    return operand, s


def reduce_sqrt2(operand, s=1):
    factors = factor(operand)

    if factors is None or factors == {operand: 1}:
        return operand, s

    for f in factors:
        if factors[f] > 1:  # if there is a square or greater in the factors, push it to the coeff
            power = pow(f, factors[f] // 2)
            s *= power
            operand //= (power * power)
    return operand, s


def test_reduce_sq():
    time_1, time_2 = 0, 0
    for _ in range(500):
        n = randrange(5, 5000)
        time_1 += timeit(lambda: reduce_sqrt1(n), number=1000)
        time_2 += timeit(lambda: reduce_sqrt2(n), number=1000)
    print(f"time for-loop1: {time_1:.2f}s")
    print(f"time for-loop2: {time_2:.2f}s")
    print(f"efficiency: {(time_1 / time_2) * 100:.2f}%")


def test_constructor_timeit():
    time_sr, time_fr = 0, 0
    n_nums = pow(10, 3)
    iters = pow(10, 3)
    for _ in range(n_nums):
        v1, v2 = randrange(2, 10), randrange(2, 10)
        time_sr += timeit(lambda: Sqrt(v1, v2, _normalize=False), number=iters)
        time_fr += timeit(lambda: fractions.Fraction(v1, v2), number=iters)

    print(f"time sr: {time_sr:.2f}s")
    print(f"time fr: {time_fr:.2f}s")
    print(f"average construction time per instance: {pow(10, 6) * (time_sr / (n_nums * iters)):.2f}µs")


if __name__ == '__main__':
    test_all()
