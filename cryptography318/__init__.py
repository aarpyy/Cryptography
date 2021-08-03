from .crypto_functions import (toBase, fromBase, StringToNum, NumToString, ApplyMult, ExtendedGCD, GCD, ModularInverse,
                               makeChineseRemainder, ChineseRemainder, baby_step_giant_step, pohlig_hellman, DSA, PollardP1, FactorInt,
                               QuadraticSieve)
from .bailliepsw_helper import Jacobi
from .prime import IsPrime, MillerRabin_bases, RandomPrime, ConfirmPrime, NextPrime, PrevPrime
from .linear_algebra import Matrix, LinearMap
import setuptools.version


__version__ = setuptools.version.__version__
__author__ = 'Andrew Carpenter'
