import math
from cryptography318.linear_algebra_deprecated import *
from cryptography318.matrix_deprecated import *
from cryptography318.linear_algebra import *
from cryptography318.prime import *
from cryptography318.crypto_functions import *
import numpy
from functools import reduce
import time
from itertools import combinations

g = 5
h = 676468030
p = 717994967
exp = 13
B = 20
# while True:
#     smooth_nums = []
#     c = pow(g, exp, p)
#     while len(smooth_nums) < 10:
#         exp += 1
#         c *= g
#         if c > p:
#             c %= p
#
#         if BSmoothQ(c, B):
#             smooth_nums.append(exp)
# smooth_nums = [2966, 2967, 22222, 22223, 51583, 51584, 51585, 51586, 51587, 102864]

smooth_nums = [2966, 2967, 22222, 22223, 51583, 51584, 51585, 51586, 51587, 102864]
mat = []
for i in smooth_nums:
    factors = factor_base_exp(pow(g, i, p), PrimesLT(B))
    factors.append(i)
    mat.append(factors)
system = Matrix(mat, aug=True)
system = system.astype(numpy.float64)


def reduce(c, r):
    for i in range(len(system)):
        row = system[i]
        if i != r and row[c] != 0:
            row -= (row[c] * system[r])

#
# system[9] /= 8
# reduce(0, 9)
# reduce(1, 8)
# reduce(2, 1)
# system[2] /= -2
# reduce(3, 2)
# system[0] /= 2
# reduce(4, 0)

print(system.rref())

logs = {2: 386603800, 3: 584671346, 5: 1, 7: 138593831, 11: 183202078, 13: 163170104, 17: 697282118, 19: 8417856}
primes = PrimesLT(B)

# k = 0
# while True:
#     k += 1
#     x = (h * pow(g, -k, p)) % p
#     if BSmoothQ(x, B):
#         exponents = factor_base_exp(x, primes)
#         # reduce sums list returned by map, which returns list of products of exponents and logs
#         print(reduce(lambda a, b: a + b, list(map(lambda i, n: i * logs[n], exponents, logs))) + k)
#         break

# assert pow(g, 2881663596, p) == h % p
