from cryptography318.core.fraction import Fraction
from cryptography318.core.expr import *
from cryptography318.core.sqrt import Sqrt
from math import sqrt


def test_constructor():
    a = Add(15, Mul(3, Sqrt(3)))
    f = Fraction(a, 3)
    assert f == Add(5, Sqrt(3))
    assert f == 5 + sqrt(3)

    f = Fraction(a, Add(5, Mul(Sqrt(3))))
    assert isinstance(f.numerator, int) and f == 3


if __name__ == '__main__':
    f = Fraction(2, Sqrt(3))
    print(f)
    print(Fraction(Sqrt(3), 4))
    print(Fraction(Sqrt(6), Sqrt(3)))
