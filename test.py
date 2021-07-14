import random

from crypto_functions import *

# p = 822171430891193443624252802027
# q = 411085715445596721812126401013
# g = 3
# a = 815739646779227
# A = pow(g, a, p)
# D = StringToNum("This is my sample document D.")
# k = random.randrange(pow(2, 10))
# S1 = pow(g, k, p) % q
# S2 = ((D + a*S1) * pow(k, -1, q)) % q
#
# print(DSA(D, S1, S2, g, p, q, A))

# for _ in range(50):
#     a, b = random.randrange(2, 100), random.randrange(2, 100)
#     if not MyExtendedGCD(a, b):
#         print(a, b)
#     elif ExtendedGCD(a, b) != MyExtendedGCD(a, b):
#         print(a, b)

# print(MyExtendedGCD(13, 78), MyExtendedGCD(78, 13))
# testGCD(39, 14)
ExtendedGCD(39, 14)
