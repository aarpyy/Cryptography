from cryptography318.utils.utils import *
import numpy as np
import pytest


@pytest.mark.parametrize('args', [
    ([1, 0, 1, 2, 0, 3], (0, 2, 3, 5)),
    ([[1, 2, 3], [0, 0], 5, [[6, 0, [0, 1]]]], ())
])
def test_where(args):
    a, b = args
    print(where(a))
    print(np.where(a))
    assert all(x == y for x, y in zip(np.where(a)[0], where(a)))


@pytest.mark.parametrize('a', [
    [[1, 2, 3], 4, [5, [6, 7]]],
    [[1, 2], [3, 4], [5, 6]],
    [[1, 2], [[3, 4], 5], [6, 7]],
    3
])
def test_shape(a):
    assert all(x == y for x, y in zip(shape(a), np.shape(a)))

