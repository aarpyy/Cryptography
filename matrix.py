import numpy


def MakeMatrix(rows, cols):
    # gets input from user for rows of matrix to be made
    # splits into imaginary and real parts when necessary
    mat = []
    for i in range(rows):
        mat.append([])
        for j in range(cols):
            n = input(f"{i + 1},{j + 1}: ")
            if "," in n:
                c = n.split(",")
                real = float(c[0])
                imag = float(c[1])
                z = complex(real, imag)
                mat[i].append(z)
            else:
                mat[i].append(float(n))
    return numpy.array(mat)


def MultiplyMatrix(matrixA=None, matrixB=None):
    if matrixA is None or matrixB is None:
        print("Enter dimensions for Matrix A: ")
        rowsA = int(input("Rows: "))
        colsA = int(input("Columns: "))

        matrixA = MakeMatrix(rowsA, colsA)

        print("Enter dimensions for Matrix B: ")
        print(f"Rows: {colsA}")
        colsB = int(input("Columns: "))

        matrixB = MakeMatrix(colsA, colsB)

    M = numpy.matmul(matrixA, matrixB)
    return M


def SquareMatrix(matrix=None):
    if matrix is None:
        print("Enter dimensions for Matrix: ")
        rows = int(input("Rows: "))
        cols = int(input("Cols: "))

        matrix = MakeMatrix(rows, cols)

    return numpy.matmul(matrix, matrix)


def InvertMatrix(matrix=None):
    if matrix is None:
        size = int(input("Enter matrix dimensions nxn: "))
        matrix = MakeMatrix(size, size)

    return numpy.linalg.inv(matrix)


def MatrixFloat(matrix):
    return matrix.astype(float)


def MatrixInt(matrix):
    return matrix.astype(int)


def ArrayToList(matrix):
    lst = []
    for i in range(len(matrix)):
        lst.append([])
        for j in range(len(matrix[0])):
            lst[i].append(matrix[i][j])
    return lst


def SortRREF(matrix):
    """Helper function for RREF that takes in RREF with unordered rows and returns
    Matrix in true RREF"""

    M = ArrayToList(matrix)
    pivot_index = 0
    for i in range(len(matrix[0]) - 1):
        index = -1
        for j in range(len(matrix)):
            if M[j][i] == 1:
                index = j
        if index != -1:
            M[index][:], M[pivot_index][:] = M[pivot_index][:], M[index][:]
            pivot_index += 1

    return numpy.array(M)


def RREF(matrix=None):
    if matrix is None:
        print("Enter dimensions for Matrix: ")
        rows = int(input("Rows: "))
        cols = int(input("Cols: "))
        matrix = MakeMatrix(rows, cols)

    row_length = len(matrix[0])
    col_length = len(matrix)
    pivot_rows = []
    for i in range(row_length):
        pivot = False
        for j in range(col_length):
            if matrix[j][i] != 0 and j not in pivot_rows:
                shift = 1 / matrix[j][i]
                for k in range(row_length):
                    matrix[j][k] *= shift
                pivot_rows.append(j)
                pivot = True

            if pivot:
                for k in range(col_length):
                    if k != j and matrix[k][i] != 0:
                        shift = matrix[k][i]
                        for l in range(row_length):
                            matrix[k][l] -= shift * matrix[j][l]
                break

    return SortRREF(matrix)
