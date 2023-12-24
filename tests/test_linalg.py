from cryptography318.linalg.linalg import *
from timeit import timeit
import numpy as np
import pytest
import random


def test_rref():
    import sympy
    np.random.seed(10)

    # Matches Sympy for random matrices
    for _ in range(200):
        a = np.random.random((50, 10)) * 10
        r = rref(a)
        assert np.allclose(np.array(sympy.Matrix(a).rref()[0], dtype=np.float64), r), print(repr(a))

    # Also test a handful that purposefully have the third column all zeros
    for _ in range(50):
        a = np.random.random((50, 10)) * 10
        a[:, 2] = np.zeros((a.shape[0],))
        r = rref(a)
        assert np.allclose(np.array(sympy.Matrix(a).rref()[0], dtype=np.float64), r), print(repr(a))

    # Raises ValueError for mishapen arrays
    with pytest.raises(Exception) as exc:
        rref([0])
    assert exc.type == ValueError

    with pytest.raises(Exception) as exc:
        rref([[[1]], [[2]]])
    assert exc.type == ValueError

    with pytest.raises(Exception) as exc:
        rref([[[1]], [[2]]])
    assert exc.type == ValueError


@pytest.mark.skip
def test_rref_timeit():
    np.random.seed(10)
    random.seed(10)

    n = 100
    it = 10
    t1 = t2 = 0
    for _ in range(it):
        a = np.random.random((random.randint(10, 50), random.randint(10, 50))) * 100
        # b = sympy.Matrix(a)
        t1 += timeit(lambda: rref(a), number=n)
        # t2 += timeit(lambda: b.rref(), number=n)
    print(f"\nTime t1: {(t1 / (it * n)) * 10e3}ms; time t2: {(t2 / (it * n)) * 10e3}ms")


def test_binary_kernel():
    np.random.seed(10)

    for _ in range(100):
        a = np.random.random((10, 5)) * 10
        a = a.astype(np.int8)
        b = kernel_gf2(a)
        assert all(np.sum(b @ np.transpose(a), axis=1) & 1 == 0)
        assert all(np.sum(a @ b.T, axis=1) & 1 == 0)


@pytest.mark.skip
def test_kernel_timeit():
    np.random.seed(10)
    random.seed(10)

    n = 100
    it = 10
    t1 = t2 = 0
    for _ in range(it):
        a = np.random.random((random.randint(10, 50), random.randint(10, 50))) * 100
        t1 += timeit(lambda: kernel(a), number=n)
        t2 += timeit(lambda: kernel_gf2(a), number=n)
    print(f"\nTime t1: {(t1 / (it * n)) * 10e3}ms; time t2: {(t2 / (it * n)) * 10e3}ms")


def test_kernel():
    np.random.seed(10)

    # Tall matrices
    for _ in range(100):
        a = np.random.random((50, 10)) * 100
        k = kernel(a)

        # If kernel is empty, nothing to check
        if len(k.shape) <= 1 or any(s == 0 for s in k.shape):
            continue
        assert all(np.isclose(np.sum(r), 0) for r in a @ k.T)

    # Wide matrices
    for _ in range(100):
        a = np.random.random((10, 50)) * 100
        k = kernel(a)

        # If kernel is empty, nothing to check
        if len(k.shape) <= 1 or any(s == 0 for s in k.shape):
            continue
        assert all(np.isclose(np.sum(r), 0) for r in a @ k.T)
