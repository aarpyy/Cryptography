import setuptools.version

from .factor import *
from .prime import *
from .utils import *
from .linalg import *
from .dlp import *

import sys
if sys.version_info < (3, 10):
    raise ImportError("cryptography318 requires Python 3.10 or higher!")
del sys

try:
    import numpy
    import sympy
except ImportError:
    raise ImportError("cryptography318 depends on both numpy and sympy!")
else:
    del numpy
    del sympy


__all__ = [
    "factor", "pollard_rho_factor", "pollard_p1",

    "siqs",

    "randprime", "prime_range", "next_prime", "prev_prime", "isprime", "miller_rabin", "baillie_psw", "confirm_prime",
    "sqrt_mod", "quadratic_residue", "quadratic_non_residue", "chinese_remainder",

    "ecm_mont", "ecm_weierstrass", "lenstra_ecm",

    "extended_gcd",

    "rref", "kernel", "binary_kernel", "minor", "det", "eigvals", "eigvec", "char_poly", "flatten",
    "matmul", "transpose",

    "baby_step_giant_step", "pollard_rho_dlp", "pohlig_hellman"
]

__version__ = setuptools.version.__version__
__author__ = 'Andrew Carpenter'
