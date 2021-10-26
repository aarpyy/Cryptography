from cryptography318.core.expr import *
from cryptography318.core.fraction import Fraction
from cryptography318.core.sqrt import Sqrt
from math import sqrt
from random import randrange
from timeit import timeit
from time import time


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


def test_random_exprs():
    start = time()
    n_exprs = pow(10, 3)
    iters = pow(10, 3)
    time_add_cons, time_mul_cons, time_add_smplfy, time_mul_smplfy = 0, 0, 0, 0
    for _ in range(n_exprs):
        rand_sq1, rand_sq2 = randrange(2, 100), randrange(2, 100)
        rand_num, rand_denom = randrange(2, 100), randrange(2, 100)
        rand_int = randrange(2, 100)
        sq1, sq2 = Sqrt(rand_sq1), Sqrt(rand_sq2)
        fr = Fraction(rand_num, rand_denom)
        time_add_cons += timeit(lambda: Add(sq1, sq2, fr, rand_int), number=iters)
        time_mul_cons += timeit(lambda: Mul(sq1, sq2, fr, rand_int), number=iters)
        m = Mul(sq1, sq2, fr, rand_int)
        time_mul_smplfy += timeit(lambda: m.simplify(), number=iters)
        a = Add(sq1, fr, rand_int, m)
        time_add_smplfy += timeit(lambda: a.simplify(), number=iters)

    print(f"Add construction time: {(time_add_cons / (n_exprs * iters)) * pow(10, 6):.2f}µs")
    print(f"Mul construction time: {(time_mul_cons / (n_exprs * iters)) * pow(10, 6):.2f}µs")
    print(f"Add simplify time: {(time_add_smplfy / (n_exprs * iters)) * pow(10, 6):.2f}µs")
    print(f"Mul simplify time: {(time_mul_smplfy / (n_exprs * iters)) * pow(10, 6):.2f}µs")
    print(f"total test time: {time() - start:.2f}s")


if __name__ == '__main__':
    test_random_exprs()


