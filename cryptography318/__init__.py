from .crypto_functions import (toBase, fromBase, StringToNum, NumToString, ApplyMult, ExtendedGCD, GCD, ModularInverse,
                               MakeChineseRemainder, ChineseRemainder, BSGS, PohligHellman, DSA)
from .matrix import MakeMatrix, MultiplyMatrix, SquareMatrix, InvertMatrix, RREF
from .bailliepsw_helper import Jacobi
from .prime import IsPrime, MillerRabin_bases, RandomPrime, ConfirmPrime, NextPrime, PrevPrime, PollardP1, FactorInt
import setuptools.version


__version__ = setuptools.version.__version__
__author__ = 'Andrew Carpenter'
