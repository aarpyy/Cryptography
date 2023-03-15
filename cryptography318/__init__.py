import sys

import setuptools.version

from cryptography318.dlp import *
from cryptography318.factor import *
from cryptography318.linalg import *
from cryptography318.prime import *
from cryptography318.utils import *

if sys.version_info < (3, 10):
    raise ImportError("cryptography318 requires Python 3.10 or higher!")
del sys

try:
    import numpy
except ImportError:
    raise ImportError("cryptography318 depends on numpy!")
else:
    del numpy

__all__ = [
    "factor", "pollard_rho_factor", "pollard_p1",

    "siqs",

    "randprime", "prime_range", "next_prime", "prev_prime", "isprime", "miller_rabin", "baillie_psw", "confirm_prime",
    "sqrt_mod", "quadratic_residue", "quadratic_non_residue", "chinese_remainder",

    "ecm_mont", "ecm_weierstrass", "lenstra_ecm",

    "extended_gcd",

    "rref", "kernel_gf2", "kernel",

    "baby_step_giant_step", "pollard_rho_dlp", "pohlig_hellman"
]

__version__ = setuptools.version.__version__
__author__ = 'Andrew Carpenter'
