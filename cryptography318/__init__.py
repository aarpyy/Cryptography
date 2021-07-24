from .crypto_functions import (toBase, fromBase, StringToNum, NumToString, ApplyMult, ExtendedGCD, GCD, ModularInverse,
                               MakeChineseRemainder, ChineseRemainder, BSGS, PohligHellman, DSA)
from .matrix import MakeMatrix, RandomMatrix, MultiplyMatrix, SquareMatrix, augmentMatrix, ResetType, separateMatrix
from .bailliepsw_helper import Jacobi
from .prime import IsPrime, MillerRabin_bases, RandomPrime, ConfirmPrime, NextPrime, PrevPrime, PollardP1, FactorInt
from .linear_algebra import (IsConsistent, Solve, IsSolvable, InvertMatrix, RREF, augRREF, IsRREF, augIsRREF,
                             MatrixEquals)
import setuptools.version


__version__ = setuptools.version.__version__
__author__ = 'Andrew Carpenter'
