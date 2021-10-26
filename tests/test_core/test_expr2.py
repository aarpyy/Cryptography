from cryptography318.core.expr import *
from cryptography318.core.sqrt import Sqrt
from cryptography318.core.fraction import Fraction
from time import time
from random import randrange


def test_simplify():
    start = time()
    n_exprs = pow(10, 5)
    a, m, k = None, None, Mul()
    try:
        for _ in range(n_exprs):
            rand_sq1, rand_sq2 = randrange(2, 100), randrange(2, 100)
            rand_num, rand_denom = randrange(2, 100), randrange(2, 100)
            rand_int = randrange(2, 100)
            sq1, sq2 = Sqrt(rand_sq1), Sqrt(rand_sq2)
            fr = Fraction(rand_num, rand_denom)

            m = Mul(sq1, sq2, fr, rand_int)
            k = m.simplify()

            # Fractions hash same as ints but don't reduce with ints, so Mul(21, Fraction(21, 1)) would produce
            # args = (21, Fraction(21, 1)) but set(args) = {21}
            assert len(k.dict[Fraction]) == len(set(k.dict[Fraction]))
            assert all(len(v) == len(set(v)) for k, v in k.dict.items() if k is not Fraction)

            a = Add(sq1, fr, rand_int, m)
            k = a.simplify()
            assert len(k.dict[Fraction]) == len(set(k.dict[Fraction]))
            assert all(len(v) == len(set(v)) for k, v in k.dict.items() if k is not Fraction)

    except AttributeError:
        print(f"non-Expr object was returned:\n{repr(m)}\n{repr(a)}\n{k}")
    except AssertionError:
        print(f"Expr failed to be simplified:\n{repr(m)}\n{repr(a)}\n{k}")

    _time = time() - start
    print(f"total test time: {_time:.2f}s")
    print(f"simplificaiton time (Add.simplify() & Mul.simplify()): {pow(10, 6)*(_time / n_exprs):.2f}µs")


def test_true_div():
    a = Mul(2, 3, Sqrt(2))
    b = Mul(Sqrt(2))
    assert a / b == Mul(6)

    a = Mul(10)
    b = Mul(2, Sqrt(5))
    assert a / b == Mul(Sqrt(5))

    a = Mul(10, Sqrt(2), Sqrt(3))
    assert a / 5 == Mul(2, Sqrt(2), Sqrt(3))

    a = Mul(10, Sqrt(2), Sqrt(3))
    assert a / Sqrt(3) == Mul(10, Sqrt(2))

    a = Mul(9, Sqrt(2), Sqrt(5), Fraction(2, 3))
    assert a / Sqrt(3) == Mul(Sqrt(2), Sqrt(3), Sqrt(5), Fraction(2, 3), 3)

    a = Mul(9, Sqrt(15), Sqrt(2))
    assert a / Sqrt(3) == Mul(9, Sqrt(2), Sqrt(5))

    a = Mul(Sqrt(3), Sqrt(6))
    assert 18 / a == Mul(Sqrt(3), Sqrt(6))
    assert Add(18) / a == Add(Mul(Sqrt(3), Sqrt(6)))

    a = Add(15, Mul(3, Sqrt(2)), Fraction(15, 2))
    b = Add(5, Fraction(5, 2), Sqrt(2))
    assert a / b == 3

    assert a / 3 == Add(5, Sqrt(2), Fraction(5, 2))


def test_gcd():
    a = Add(5, Mul(10, Sqrt(2)), Mul(15, Sqrt(3)))
    assert a.gcd == 5

    a += 5
    assert a.gcd == 5

    a += Mul(5, Sqrt(3))
    assert a.gcd == 10

    b = Mul(13, Sqrt(2), Sqrt(3))
    assert b.gcd == 13


if __name__ == '__main__':
    print(Mul(Sqrt(3)) * Sqrt(6))
    print(Sqrt(3) * Sqrt(6))
