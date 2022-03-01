import setuptools.version

from factor import factor, pollard_rho_factor, pollard_p1
from siqs import siqs
from prime import randprime, prime_range, next_prime, prev_prime
from elliptic import ecm_mont, ecm_weierstrass, lenstra_ecm
from utils import extended_gcd
from linalg import (rref, kernel, binary_kernel, minor, det, eigvals, eigvec, char_poly, flatten,
                    matmul, vec_matmul, transpose)
from dlp import baby_step_giant_step, pollard_rho_dlp, pohlig_hellman

__all__ = [
    "factor", "pollard_rho_factor", "pollard_p1",

    "siqs",

    "randprime", "prime_range", "next_prime", "prev_prime",

    "ecm_mont", "ecm_weierstrass", "lenstra_ecm",

    "extended_gcd",

    "rref", "kernel", "binary_kernel", "minor", "det", "eigvals", "eigvec", "char_poly", "flatten",
    "matmul", "vec_matmul", "transpose",

    "baby_step_giant_step", "pollard_rho_dlp", "pohlig_hellman"

]

__version__ = setuptools.version.__version__
__author__ = 'Andrew Carpenter'
