import setuptools.version

from .factor import factor, pollard_rho_factor, pollard_p1
from .siqs import siqs
from .prime import randprime, prime_range, next_prime, prev_prime


__version__ = setuptools.version.__version__
__author__ = 'Andrew Carpenter'
