from cryptography318.linear_algebra import *
from cryptography318.crypto_functions import *



m = numpy.array([[4, 2, 2],
                 [4, 6, 6]])

r = m.copy().astype(numpy.float128)
print(isinstance(r[0][0], float))
print(type(r[0][0]))
