import sys

import setuptools.version

from cryptography318.dlp import baby_step_giant_step, pollard_rho_dlp, pohlig_hellman
from cryptography318.factor import pollard_rho_factor, pollard_pm1, qs, factor, ecm
from cryptography318.linalg import kernel_gf2
from cryptography318.prime import (baillie_psw, chinese_remainder, confirm_prime, isprime, lift_sqrt, miller_rabin,
                                   next_prime, prev_prime, quadratic_non_residue, quadratic_residue,
                                   randprime, sqrt_mod, primesieve)
from cryptography318.utils import extended_gcd

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
    "factor", "pollard_rho_factor", "pollard_pm1",

    "qs",

    "randprime", "next_prime", "prev_prime", "isprime", "miller_rabin", "baillie_psw", "confirm_prime",
    "sqrt_mod", "quadratic_residue", "quadratic_non_residue", "chinese_remainder", "lift_sqrt", "primesieve",

    "ecm",

    "extended_gcd",

    "kernel_gf2",

    "baby_step_giant_step", "pollard_rho_dlp", "pohlig_hellman"
]

__version__ = setuptools.version.__version__
__author__ = 'Andrew Carpenter'
