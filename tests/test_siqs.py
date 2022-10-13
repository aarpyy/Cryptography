from cryptography318.prime import randprime
from cryptography318.factor.siqs import siqs
from cryptography318.factor.siqs_np import siqs as siqs_np
from cryptography318.factor.siqs_np import smooth_t, smooth_u, primes, solve_matrix
from cryptography318.linalg import binary_kernel_np
import numpy as np

from pathlib import Path


test_root = Path(__file__).parent.absolute()


def test_solve_matrix():
    global smooth_t, smooth_u, primes
    smooth_u = []
    with open(test_root.parent.joinpath("smooth_u.txt"), "r") as fp:
        for line in fp.readlines():
            smooth_u.append([int(v) for v in line.strip('\n').split(",")])

    smooth_t = []
    with open(test_root.parent.joinpath("smooth_t.txt"), "r") as fp:
        for line in fp.readlines():
            smooth_t.append([int(v) for v in line.strip('\n').split(",")])

    primes = []
    with open(test_root.parent.joinpath("primes.txt"), "r") as fp:
        for line in fp.readlines():
            primes.append([int(v) for v in line.strip('\n').split(",")])

    n = 26910013076074612954386739
    solve_matrix(n)


def test_siqs_34():
    lower = pow(10, 12)
    upper = lower * 10
    a = randprime(lower, upper)
    b = randprime(lower, upper)
    print(f"Factors: {{a: {a}, b: {b}}}")
    c = siqs_np(a * b, fp=str(test_root.parent.joinpath("cryptography318").joinpath("prime/primes.txt")), loud=True)
    assert c == a or c == b


if __name__ == "__main__":
    test_siqs_34()
