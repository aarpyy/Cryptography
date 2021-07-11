from crypto_functions import *

for _ in range(500):
    n = random.randrange(2, 10)
    A = []
    for row in range(n):
        A.append([])
        for col in range(n):
            A[row].append(random.randrange(-10, 10))
    B = []
    for row in range(n):
        B.append([])
        for col in range(n):
            B[row].append(random.randrange(-10, 10))

    A = numpy.array(A)
    B = numpy.array(B)

    A2B2 = MultiplyMatrix(SquareMatrix(A), SquareMatrix(B))
    AB2 = SquareMatrix(MultiplyMatrix(A, B))
    count = 0
    for i in range(len(A2B2)):
        for j in range(len(A2B2[i])):
            if A2B2[i][j] != AB2[i][j]:
                count += 1
                break
    if count == 0:
        print(A)
        print(B)
