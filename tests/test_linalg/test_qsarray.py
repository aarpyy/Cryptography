from cryptography318.linalg.qsarray import *
from random import randrange
from cryptography318.linalg.linalg import BinaryMatrix, ArrayMod
from timeit import timeit
from cryptography318.numbers.quadratic_sieve import *
from cryptography318.numbers.prime import *

import numpy as np


def test_add():
    a = barray(randrange(32), _len=5)
    b = barray(randrange(32), _len=5)
    c = a + b
    for i, bit in enumerate(c):
        assert (a[i] + b[i]) & 1 == bit


def test_mul():
    a = barray(13, _len=5)
    b = barray(21, _len=5)
    # dot product of 01101 and 10101 = 1*1 + 1*1
    assert a * b == 2

    a = barray(17, _len=5)
    b = barray(30, _len=5)
    assert a * b == 1


def test_get_item():
    a = barray(21, _len=5)
    print(a)
    print(a[2:])
    print(a[:2])


def test_matrix():
    a = bmatrix([barray(31, _len=5), barray(17, _len=5), barray(8, _len=5)])
    print(a)
    # b = a.kernel()
    a.kernel()


def test_kernel():
    a = bmatrix([barray(31, _len=5), barray(17, _len=5), barray(8, _len=5), barray(16, _len=5), barray(8, _len=5),
                 barray(4, _len=5), barray(2, _len=5), barray(1, _len=5)]).transpose()

    _iters = pow(10, 4)

    time_b = timeit(lambda: a.transpose(), number=_iters)
    as_list = list(map(lambda r: ArrayMod(r, 2), list(map(list, a))))
    m = BinaryMatrix(as_list)
    n = np.array(as_list)
    time_bm = timeit(lambda: m.transpose(), number=_iters)
    time_np = timeit(lambda: n.transpose(), number=_iters)
    time_py = timeit(lambda: list(map(list, zip(*as_list))), number=_iters)
    print(f"time bmatrix: {time_b:.3f}")
    print(f"time BinaryMatrix: {time_bm:.3f}")
    print(f"time numpy: {time_np:.3f}")
    print(f"time python: {time_py:.3f}")


def test_all():
    test_add()
    test_mul()


def separate(arr, j):
    left, right = [], []
    for i in range(len(arr)):
        left.append(arr[i][:j])
        right.append(arr[i][j:])
    return left, right


if __name__ == '__main__':
    N = 3023483
    print(quadratic_sieve(N))
