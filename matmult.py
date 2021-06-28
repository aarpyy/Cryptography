from crypto_functions import *

A = InvertMatrix()
B = MakeMatrix(3, 1)
M = numpy.matmul(A, B)
for row in M:
    for i in range(len(row)):
        row[i] = row[i] * 5
print(M)



