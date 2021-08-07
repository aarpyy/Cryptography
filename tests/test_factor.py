from cryptography318.crypto_functions import _factorPerfectSquare, SolveDLP
from cryptography318.prime import *
from cryptography318.crypto_functions import *
from math import prod
import pytest


@pytest.mark.skip
def test_factor_perfect_square():
    factors = _factorPerfectSquare(198103, B=20)
    n = 1
    for f in factors:
        n *= pow(f, factors[f])
    assert n == 198103


@pytest.mark.skip
def test_factor_int(power=1):
    from itertools import combinations as _all

    def same_combination(lst1, lst2):
        not_same = 0
        for term in lst1:
            in_other = False
            for term1 in lst2:
                same = 0
                for t in term:
                    if t not in term1:
                        continue
                    same += 1
                if same == len(term):
                    in_other = True
                    break
            if not in_other:
                not_same += 1
        return not_same == 0
    import timeit
    p = RandomPrime(pow(10, power), pow(10, power + 1)) * RandomPrime(pow(10, power + 1), pow(10, power + 2))
    # factors = FactorInt(p)
    # assert p == prod(map(lambda b, e: pow(b, e), factors.keys(), factors.values()))
    print(timeit.timeit(lambda: QuadraticSieve(p), number=100))
    print(timeit.timeit(lambda: FactorInt(p), number=100))


def test_pollard_rho(it=50):
    for _ in range(it):
        g = 4
        p = RandomPrime(pow(2, 30))
        e = randrange(p - 1)
        h = pow(g, e, p)
        r = SolveDLP(g, h, p)
        assert pow(g, r, p) == h


@pytest.mark.skip
def test_all_tests():
    test_factor_perfect_square()
    test_factor_int()


if __name__ == '__main__':
    test_pollard_rho()
