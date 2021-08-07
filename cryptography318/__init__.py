from .crypto_functions import (toBase, fromBase, StringToNum, NumToString, ExtendedGCD, makeChineseRemainder,
                               ChineseRemainder, baby_step_giant_step, pohlig_hellman, DSA,
                               pollard_p1, FactorInt, index_calculus_dlp, StringToElliptic, EllipticToString, Elliptic,
                               EllipticCurve, SolveDLP, elliptic_bsgs, calculate_state)
from .bailliepsw_helper import Jacobi
from .prime import IsPrime, MillerRabin_bases, RandomPrime, ConfirmPrime, NextPrime, PrevPrime, BailliePSW_Primality
from .linear_algebra import Matrix, LinearMap, where, aslist, isnumber, all_elements
import setuptools.version


__version__ = setuptools.version.__version__
__author__ = 'Andrew Carpenter'
