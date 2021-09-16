from cryptography318.numbers.quadratic_sieve import *
from timeit import timeit


def test_efficiency():
    times = {}
    curr_time = 0
    try:
        for i in range(pow(2, 10), pow(2, 15)):
            curr_time += timeit(lambda: quadratic_sieve(i, force=20), number=10)
            if i % 10000 == 0:
                times[i] = curr_time
                curr_time = 0
        print(times)
    except IndexError:
        print(i, times)


if __name__ == '__main__':
    test_efficiency()
