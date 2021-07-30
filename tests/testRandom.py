from cryptography318.linear_algebra_deprecated import *
from cryptography318.matrix_deprecated import *
from cryptography318.linear_algebra import *
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

print(T.change_basis(b))



# print("3.5.8 a)")
# print(f"i) \n{ChangeBasis(a1, b)}\n")
# print(f"ii) \n{ChangeBasis(a2, b)}\n")
# print(f"iii) \n{ChangeBasis(a3, c)}\n")
# print(f"iv) \n{ChangeBasis(a4, c)}\n")
# print(f"v) \n{ChangeBasisMap(T, in_basis=c)}\n")
# print(f"vi) \n{ChangeBasisMap(T, out_basis=b)}\n")
# print(f"vii) \n{ChangeBasisMap(T, in_basis=c, out_basis=b)}\n")
