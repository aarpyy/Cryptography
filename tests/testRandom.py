import math

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

# print(T.change_basis(b))

m = Matrix([[complex(0, 1), complex(4, -2)],
            [complex(1, 2), complex(2, 1)]])
print(m)
print(T)



# print("3.5.8 a)")
# print(f"i) \n{ChangeBasis(a1, b)}\n")
# print(f"ii) \n{ChangeBasis(a2, b)}\n")
# print(f"iii) \n{ChangeBasis(a3, c)}\n")
# print(f"iv) \n{ChangeBasis(a4, c)}\n")
# print(f"v) \n{ChangeBasisMap(T, in_basis=c)}\n")
# print(f"vi) \n{ChangeBasisMap(T, out_basis=b)}\n")
# print(f"vii) \n{ChangeBasisMap(T, in_basis=c, out_basis=b)}\n")




# imag = n.imag
#                     if n.imag < 1:
#                         imag = n.imag * -1
#                         str_imag += "-"
#                     else:
#                         str_imag += "+"
#                     if n.real == 0 and imag == 0:
#                         str_real = "0"
#                     elif n.real == 0 and (imag == 1 or imag == -1):
#                         str_imag += " i"
#                     elif n.real == 0:
#                         s = str(imag).split(".")
#                         if s[1] == '0':
#                             str_imag += s[0] + "i"
#                         else:
#                             str_imag += str(imag) + "i"
#                         str_imag = str_imag[2:]
#                     elif n.imag == 0:
#                         s = str(n.real).split(".")
#                         if s[1] == '0':
#                             str_real = s[0]
#                         else:
#                             str_real = str(n.real)
#                     else:
#                         s = str(n.real).split(".")
#                         if s[1] == '0':
#                             str_real = s[0]
#                         else:
#                             str_real = str(n.real)
#                         s = str(imag).split(".")
#                         if imag == 1:
#                             str_imag += " i"
#                         elif s[1] == '0':
#                             str_imag += " " + s[0] + "i"
#                         else:
#                             str_imag += " " + str(imag) + "i"
#                     string = str_real + " " + str_imag
#                     if len(string) > max_len:
#                         max_len = len(string)
#                     str_array[i].append(string)
