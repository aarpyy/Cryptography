from cryptography318.linalg.linalg import *
import numpy as np
import random


def test_binary_kernel():
    np.random.seed(10)

    for _ in range(100):
        a = np.random.random((10, 5)) * 10
        a = a.astype(np.int8)
        b = kernel_gf2(a)
        assert all(np.sum(b @ np.transpose(a), axis=1) & 1 == 0)
        assert all(np.sum(a @ b.T, axis=1) & 1 == 0)
