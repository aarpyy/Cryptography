from cryptography318.core.fraction import *
from fractions import Fraction as PyFraction
from timeit import timeit
from random import randrange
from numbers import *
from math import sqrt
from cryptography318.core.sqrt import Sqrt
import pytest


def test_constructor_time():
    n_nums, n_iter = pow(10, 3), pow(10, 3)
    time_fr, time_pyf = 0, 0
    for _ in range(n_nums):
        a, b = randrange(2, 100), randrange(2, 100)
        time_fr += timeit(lambda: Fraction(a, b), number=n_iter)
        time_pyf += timeit(lambda: PyFraction(a, b), number=n_iter)

    print(f"fraction took: {time_fr:.2f}s")
    print(f"pyfraction took: {time_pyf:.2f}s")
    print(f"average construction time per instance: {pow(10, 6) * (time_fr / (n_nums * n_iter)):.2f}µs")


def test_constructor():
    x = Fraction(2, 5)
    assert isinstance(x, Fraction) and isinstance(x, Rational)
    assert isinstance(x.numerator, int) and isinstance(x.denominator, int)
    assert x == .4
    assert x == PyFraction(2, 5)

    x = Fraction(8, 6)
    assert isinstance(x, Fraction) and isinstance(x, Rational)
    assert isinstance(x.numerator, int) and isinstance(x.denominator, int)
    assert x.numerator == 4 and x.denominator == 3

    x = Fraction(2, 1)
    assert isinstance(x, Fraction) and isinstance(x, Rational)
    assert isinstance(x.numerator, int) and isinstance(x.denominator, int)
    assert x == 2 and x.denominator == 1

    x = Fraction(2.5)
    assert isinstance(x, Fraction) and isinstance(x, Rational)
    assert isinstance(x.numerator, int) and isinstance(x.denominator, int)
    assert x == 2.5 and x.numerator == 5 and x.denominator == 2

    x = Fraction(Sqrt(2), 2)
    assert isinstance(x, float) and abs(x - sqrt(2)/2) < pow(10, -6)

    x = Fraction(Sqrt(2), 3)
    assert isinstance(x, Fraction) and isinstance(x, Rational)
    assert isinstance(x.numerator, Sqrt) and isinstance(x.denominator, int)
    assert x == sqrt(2) / 3

    with pytest.raises(TypeError) as exc:
        x = Fraction(complex(1, 0))
    assert "must be rational" in str(exc.value)

    x = Fraction('2/5')
    assert isinstance(x, Fraction) and isinstance(x, Rational)
    assert isinstance(x.numerator, int) and isinstance(x.denominator, int)
    assert x == .4

    x = Fraction('2 / 5')
    assert isinstance(x, Fraction) and isinstance(x, Rational)
    assert isinstance(x.numerator, int) and isinstance(x.denominator, int)
    assert x == .4

    x = Fraction('sqrt(2) / 5')
    assert isinstance(x, Fraction) and isinstance(x, Rational)
    assert isinstance(x.numerator, Sqrt) and isinstance(x.denominator, int)
    assert x == sqrt(2) / 5

    x = Fraction('1.25E2')
    assert isinstance(x, Fraction) and isinstance(x, Rational)
    assert isinstance(x.numerator, int) and isinstance(x.denominator, int)
    assert x == 125

    x = Fraction('-1.25E-1')
    assert isinstance(x, Fraction) and isinstance(x, Rational)
    assert isinstance(x.numerator, int) and isinstance(x.denominator, int)
    assert x == -1/8


def test_add():
    a = Fraction(2, 5)

    b = Fraction(3, 5)
    c = a + b
    assert c == 1 and isinstance(c, Fraction)

    b = Fraction(Sqrt(2, 2), 5)
    c = a + b
    assert abs(c - (2 + 2 * sqrt(2)) / 5) < pow(10, -6) and isinstance(c, float)

    c = a + 2
    assert isinstance(c, Fraction) and c.numerator == 12 and c.denominator == 5

    c = a + 2.0
    assert isinstance(c, Fraction) and c.numerator == 12 and c.denominator == 5

    # while this does result in an argument that could be a fraction, adding with non-integer floats should
    # never return a new fraction
    c = a + 2.1
    assert isinstance(c, float) and c == 2.5

    a = Fraction(Sqrt(5), 7)

    b = Fraction(Sqrt(5), 9)
    c = a + b
    assert isinstance(c, Fraction) and isinstance(c.numerator, Sqrt) and c.denominator == 63


def test_sub():
    a = Fraction(2, 5)

    b = Fraction(1, 5)
    c = a - b
    assert c == .2 and isinstance(c, Fraction)

    b = Fraction(Sqrt(2, 2), 5)
    c = a - b
    assert abs(c - (2 - 2 * sqrt(2)) / 5) < pow(10, -6) and isinstance(c, float)

    c = a - 2
    assert isinstance(c, Fraction) and c.numerator == -8 and c.denominator == 5

    c = a - 2.0
    assert isinstance(c, Fraction) and c.numerator == -8 and c.denominator == 5

    # while this does result in an argument that could be a fraction, adding with non-integer floats should
    # never return a new fraction
    c = a - 1.4
    assert isinstance(c, float) and c + 1.0 < pow(10, -6)

    a = Fraction(Sqrt(5, 2), 7)

    b = Fraction(Sqrt(5, 7), 9)
    c = a - b
    assert isinstance(c, Fraction) and isinstance(c.numerator, Sqrt) and c.denominator == 63


def test_mul():
    a = Fraction(2, 5)

    b = Fraction(3, 7)
    c = a * b
    assert isinstance(c, Fraction) and c.numerator == 6 and c.denominator == 35

    c = a * 2
    assert isinstance(c, Fraction) and c.numerator == 4 and c.denominator == 5

    c = a * 2.0
    assert isinstance(c, Fraction) and c.numerator == 4 and c.denominator == 5

    b = Sqrt(2, 2)
    c = a * b
    assert isinstance(c, Fraction) and c.numerator == Sqrt(2, 4) and c.denominator == 5

    a = Fraction(Sqrt(2, 2), 5)

    b = Fraction(2, 5)
    c = a * b
    assert isinstance(c, Fraction) and c.numerator == Sqrt(2, 4) and c.denominator == 25

    b = Fraction(Sqrt(2, 5), 5)
    c = a * b
    assert isinstance(c, Fraction) and isinstance(c.numerator, int)
    assert c.numerator == 4 and c.denominator == 5

    b = Fraction(Sqrt(2, 7), 5)
    c = a * b
    assert isinstance(c, Fraction) and isinstance(c.numerator, int)
    assert c.numerator == 28 and c.denominator == 25


def test_div():
    a = Fraction(2, 5)

    b = Fraction(5, 2)
    c = a / b
    assert isinstance(c, Fraction) and c.numerator == 4 and c.denominator == 25

    b = Fraction(2, 5)
    c = a / b
    assert isinstance(c, Fraction) and c.numerator == 1 and c.denominator == 1

    c = a // b
    assert isinstance(c, int) and c == 1

    b = Fraction(Sqrt(17))
    c = 1 / b
    assert isinstance(c, Fraction) and isinstance(c.numerator, Sqrt)
    assert c.numerator == Sqrt(17) and c.denominator == 17


def test_mod():
    a = Fraction(7, 5)

    b = Fraction(2, 5)
    c = a % b
    assert isinstance(c, Fraction) and c.numerator == 1 and c.denominator == 5

    # modding by square root should just return float even in fraction
    b = Fraction(Sqrt(2, 2), 5)
    c = a % b
    assert isinstance(c, float)


@pytest.mark.skip
def test_all():
    test_constructor()
    test_add()
    test_sub()
    test_mul()
    test_div()
    test_mod()


if __name__ == '__main__':
    test_all()
