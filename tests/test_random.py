import math

from cryptography318.linear_algebra_deprecated import *
from cryptography318.matrix_deprecated import *
from cryptography318.linear_algebra import *
from cryptography318.prime import *
from cryptography318.crypto_functions import *
import numpy

b = Matrix([[1, 1],
            [1, -1]])
c = Matrix([[1, 1, 0],
            [1, 0, 1],
            [0, 1, 1]])
a1 = Matrix([[1],
             [2]])
a2 = Matrix([[2],
             [3]])
a3 = Matrix([[1],
             [-2],
             [3]])
a4 = Matrix([[3],
             [2],
             [1]])

T = Matrix([[0, 1, 0],
            [1, 0, 0]])

# print(T.change_basis(b))


p = RandomPrime(pow(2, 25))
print(p)
print(test_cprime(p))

