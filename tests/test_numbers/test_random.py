from cryptography318.core.fraction import Fraction
from cryptography318.core.expr import *


if __name__ == '__main__':
    a = Mul(Sqrt(2), 5, 7)
    b = Add(2, Sqrt(2), Fraction(3, 4))
    c = Mul(Sqrt(8), Sqrt(2), 7, Fraction(5, 4), Fraction(7, 8), Sqrt(25), Sqrt(50))
    d = a * a
    print(d)
    d.simplify(_update=True)
    print(d, type(d))
