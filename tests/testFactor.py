from cryptography318.crypto_functions import _factorPerfectSquare
from itertools import combinations_with_replacement

# x = list(map(list, list(combinations_with_replacement([1, 2, 3, 4], 2))))

x = _factorPerfectSquare(198103, B=20)
print(x)
