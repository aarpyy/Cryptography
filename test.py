import numpy

def MakeMatrix(rows, cols):
    mat = []
    for i in range(rows):
        mat.append([])
        for j in range(cols):
            n = input(f"{i+1},{j+1}: ")
            if "," in n:
                c = n.split(",")
                real = int(c[0])
                imag = int(c[1])
                z = complex(real, imag)
                mat[i].append(z)
            else:
                mat[i].append(int(n))
    return numpy.array(mat)

print("Enter dimensions for Matrix A:")
rowsA = int(input("Rows: "))
colsA = int(input("Columns: "))

A = MakeMatrix(rowsA, colsA)

print("Enter dimensions for Matrix B: ")
print(f"Rows: {colsA}")
colsB = int(input("Columns: "))

B = MakeMatrix(colsA, colsB)

M = numpy.matmul(A, B)
print(M)
