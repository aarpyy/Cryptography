from cryptography318.numbers.prime import *
from cryptography318.numbers.crypto_functions import *
from cryptography318.numbers.quadratic_sieve import *
from cryptography318.numbers.factor import *
import time
import timeit
from timeit import timeit
import pytest


@pytest.mark.skip
def test_factor_int(power=1):
    p = randprime(pow(10, power), pow(10, power + 1)) * randprime(pow(10, power + 1), pow(10, power + 2))
    print(timeit.timeit(lambda: quadratic_sieve(p), number=100000))
    print(timeit.timeit(lambda: factor(p), number=100000))


def test_pollard_rho(it=50):
    for _ in range(it):
        g = 4
        p = randprime(pow(2, 30))
        e = randrange(p - 1)
        h = pow(g, e, p)
        r = pollard_rho_dlp(g, h, p)
        assert pow(g, r, p) == h


def test_factor_if_smooth():
    n = 5 * 7 * 7 * pow(3, 4)
    assert factor_if_smooth(n, [2, 3, 5, 7]) == [0, 4, 1, 2]
    assert factor_if_smooth(n, [2, 3, 5]) is None


def test_qs(it=10):
    for _ in range(it):
        a = randprime(pow(10, 4), pow(10, 5))
        b = randprime(pow(10, 4), pow(10, 5))
        n = a * b
        factors = quadratic_sieve(n)
        if factors is None:
            print(a, b, n)
        result = 1
        for f in factors:
            result *= pow(f, factors[f])
        assert result == n


@pytest.mark.skip
def test_max_qs():

    # this test shouldn't be run with pytest, it just purely for finding largest integer than can be factored
    # in a short amount of time by quadratic sieve
    start = time.time()
    a = randprime(pow(10, 9))
    b = randprime(pow(10, 5))
    n = a * pow(b, 2)
    print(a, b, n)
    print(quadratic_sieve(n))
    print(f"This took {time.time() - start:.2f}s")


def test_b_smooth():
    for _ in range(50):
        B = 25
        primes = primes(B)
        n = 1
        for i in range(len(primes)):
            n *= pow(primes[i], randrange(2, 5))

        assert b_smooth(n, B)
        assert b_smooth(n, factors=primes)


def jd(*args: dict) -> dict:
    """Joins multiple dictionaries in a way that sums values of shared keys. Assumes all values
    support + method."""

    def join(dict1, dict2):
        for key in dict2:
            if key in dict1:
                dict1[key] += dict2[key]
            else:
                dict1[key] = dict2[key]
        return dict1

    def update(dict1, dict2):
        dict1.update(dict2)
        return dict1

    return reduce(lambda a, b: join(a, b) if any(k in b for k in a) else update(a, b), args)


def rd_fact(factors):
    result = {}
    for f in factors:
        if not isprime(f):
            k = factor(f)
            for a in k:
                k[a] *= factors[f]
            result = join_dict(result, k)
        else:
            result = join_dict(result, {f: factors[f]})
    return result


def rd_fact2(factors):
    new_factors = []
    to_del = set()
    for f in factors:
        if not isprime(f):
            k = factor(f)
            for a in k:
                k[a] *= factors[f]
            new_factors.append(k)
            to_del.add(f)
    for f in to_del:
        del factors[f]
    for d in new_factors:
        for f in d:
            if f in factors:
                factors[f] += d[f]
            else:
                factors[f] = d[f]
    return factors


def test_compare_reduction():
    factors = {2: 4, 4: 6, 22: 7, 625: 18}
    res = {2: 23, 11: 7, 5: 72}
    iters = pow(10, 4)
    itr1 = timeit(lambda: rd_fact2(factors), number=iters)
    itr2 = timeit(lambda: rd_fact(factors), number=iters)
    print(f"for w/ new join: {itr1 * pow(10, 3):.2f}ms")
    print(f"for w/ old join: {itr2 * pow(10, 3):.2f}ms")


if __name__ == '__main__':
    n = (randprime(pow(2, 10)) + 4) * 2
    print(factor(n))
