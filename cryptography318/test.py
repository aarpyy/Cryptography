from matrix import *
from prime import *
from crypto_functions import *
from linear_algebra import *
import numpy

# m1 = RandomMatrix(rows=5, cols=3)
# m2 = RandomMatrix(rows=len(m1[0]), cols=1)
# sol = MultiplyMatrix(m1, m2)
# print(m1)
# print(RREF(m1))
# print(m2)
# print(sol)
# m3 = augmentMatrix(m1, sol)
# g = augRREF(m3)
# print(g, augIsRREF(g))
# print(Solve(m3, aug=True))

aug = numpy.array([[1., 0., 0., -46.96875],
                   [0., 1., 0., -24.],
                   [0., 0., 1., 24.],
                   [0., 0., 0., 0.],
                   [0., 0., 0., 0.]])
print(aug)
print(Solve(aug, aug=True))
