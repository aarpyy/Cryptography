from cryptography318.numbers.prime import *

range_tested = False
iter_tested = False


def test_range():
    global primesieve

    range_tested = True

    primes = []
    for p in primesieve.range(2, 20):
        primes.append(p)
    assert primes == [2, 3, 5, 7, 11, 13, 17, 19]


def test_iter():
    global primesieve
    global range_tested

    iter_tested = True

    if not range_tested:
        test_range()

    primes = []
    for p in primesieve:
        primes.append(p)

    # should still be extended past 13 by the above function
    assert primes == [2, 3, 5, 7, 11, 13, 17, 19]

    primesieve.extend(50)
    primes = []
    for p in primesieve:
        primes.append(p)
    assert primes == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]


def test_search():
    global primesieve
    global iter_tested

    if not iter_tested:
        test_iter()

    assert primesieve.search(2) == 0
    assert primesieve.search(47) == len(primesieve) - 1
    assert primesieve.search(4) == (1, 2)


if __name__ == '__main__':
    test_search()
