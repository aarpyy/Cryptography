from cryptography318.core.expr import *
from cryptography318.core.fraction import Fraction


def test_mul():
    # for all str tests, direct string equivalence is not tested since the order of Real instances
    # is not deterministic, but based on whatever order they were constructed in and could change

    a = Mul(2, 5, Sqrt(8), Fraction(5, 6))
    assert isinstance(a, Mul)
    assert 'sqrt(2)' in str(a) and '5/6' in str(a) and '20' in str(a)

    b = a * Sqrt(2)
    assert '40' in str(b) and 'sqrt(2)' not in str(b)

    b = a * Sqrt(4)
    assert 'sqrt(2)' in str(b)

    b = a * Sqrt(8)
    assert 'sqrt(2)' not in str(b)

    b = a * Fraction(6, 5)
    assert '/' not in str(b)
    # Fraction should reduce to Fraction w/ denominator == 1 which Mul should reduce to
    # int which should combine w/ _int leaving 20 * sqrt(2)
    assert str(b).count('*') == 1


if __name__ == '__main__':
    a = Mul(2, 5, Sqrt(6), Sqrt(2))
    b = Add(Mul(Sqrt(2), Sqrt(6)))
    print(a, b)
    print(a + b)

