import numpy
import random


def MakeMatrix(rows, cols, rand=False):
    """
    Function that creates and returns a numpy.array object.

    Creates matrix rows by getting input from user (or random if specified) and then appends rows to larger list.
    Complex input is handled using Python's built in complex() function. Matrix holds integers by default, only
    using floats if user enters number with decimal point.

    If rand == False, function will allow user to enter input for values of matrix. If rand == True, random
    values x will be used where -50 <= x <= 50 is a random value generated using random.randrange().
    """

    mat = []
    for i in range(rows):
        mat.append([])
        for j in range(cols):
            n = input(f"{i + 1},{j + 1}: ") if not rand else str(random.randrange(-50, 51))
            if "," in n:
                c = n.split(",")
                real = float(c[0])
                imag = float(c[1])
                z = complex(real, imag)
                mat[i].append(z)
            else:
                x = float(n) if "." in n else int(n)
                mat[i].append(x)
    return numpy.array(mat)


def RandomMatrix(rows=None, cols=None):
    if rows is None:
        rows = random.randrange(1, 10)
    if cols is None:
        cols = random.randrange(1, 10)

    return MakeMatrix(rows, cols, rand=True)


def CopyMatrix(matrix):
    new_mat = []
    for i in range(len(matrix)):
        new_mat.append([])
        for j in range(len(matrix[0])):
            new_mat[i].append(matrix[i][j])

    return numpy.array(new_mat)


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
    return numpy.linalg.inv(CopyMatrix(matrix))


def MatrixFloat(matrix, bits=16):
    float_types = {16: numpy.float16, 32: numpy.float32, 64: numpy.float64}
    if bits in float_types:
        return matrix.astype(float_types[bits])
    return matrix.astype(float)


def ResetType(matrix):
    matrix = ArrayToList(matrix)
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            e = str(matrix[i][j]).split(".")
            if len(e) > 1 and e[1] == '0':
                matrix[i][j] = int(matrix[i][j])
    return numpy.array(matrix)


def augmentMatrix(matrix, solutions):
    m = ArrayToList(matrix)
    for i in range(len(m)):
        m[i].append(solutions[i][0])
    return numpy.array(m)


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
        if index > pivot_index:
            M[index][:], M[pivot_index][:] = M[pivot_index][:], M[index][:]
        if index != -1:
            pivot_index += 1

    return numpy.array(M)


def RREF(matrix=None):
    if matrix is None:
        print("Enter dimensions for Matrix: ")
        rows = int(input("Rows: "))
        cols = int(input("Cols: "))
        matrix = MakeMatrix(rows, cols)

    matrixM = MatrixFloat(CopyMatrix(matrix))

    row_length = len(matrixM[0])
    col_length = len(matrixM)
    pivot_rows = []
    for i in range(row_length):
        pivot = False
        for j in range(col_length):
            if matrixM[j][i] != 0 and j not in pivot_rows:
                if matrixM[j][i] != 1:
                    shift = 1 / matrixM[j][i]
                    for k in range(row_length):
                        matrixM[j][k] *= shift

                pivot_rows.append(j)
                pivot = True

            if pivot:
                for k in range(col_length):
                    if k != j and matrixM[k][i] != 0:
                        shift = matrixM[k][i]
                        for l in range(row_length):
                            matrixM[k][l] -= shift * matrixM[j][l]
                break

    return ResetType(SortRREF(matrixM))


def augRREF(matrix, solutions=None):
    if solutions is not None:
        matrix = augmentMatrix(matrix, solutions)
    matrixM = MatrixFloat(CopyMatrix(matrix))

    row_length = len(matrixM[0])
    col_length = len(matrixM)
    pivot_rows = []
    for i in range(row_length - 1):
        pivot = False
        for j in range(col_length):
            if matrixM[j][i] != 0 and j not in pivot_rows:
                if matrixM[j][i] != 1:
                    shift = 1 / matrixM[j][i]
                    for k in range(row_length):
                        matrixM[j][k] *= shift

                pivot_rows.append(j)
                pivot = True

            if pivot:
                for k in range(col_length):
                    if k != j and matrixM[k][i] != 0:
                        shift = matrixM[k][i]
                        for l in range(row_length):
                            matrixM[k][l] -= shift * matrixM[j][l]
                break

    return ResetType(SortRREF(matrixM))


def IsRREF(matrix):
    """Function that returns True if a matrix is in reduced-row echelon form, and False otherwise"""

    pivots = []
    for i in range(len(matrix)):
        j = 0
        checked = False
        while j < len(matrix[0]) and not checked:
            e = matrix[i][j]
            if e not in [0, 1]:
                return False
            if e == 1:
                if j in pivots:
                    return False
                pivots.append(j)
                checked = True
            j += 1

    for c in range(len(pivots) - 1):
        if pivots[c] > pivots[c+1]:
            return False
    return True


def augIsRREF(matrix):
    """Function that returns True if a matrix is in reduced-row echelon form, and False otherwise"""

    pivots = []
    for i in range(len(matrix) - 1):
        j = 0
        checked = False
        while j < len(matrix[0]) and not checked:
            e = matrix[i][j]
            if e not in [0, 1]:
                return False
            if e == 1:
                if j in pivots:
                    return False
                pivots.append(j)
                checked = True
            j += 1

    for c in range(len(pivots) - 1):
        if pivots[c] > pivots[c+1]:
            return False
    return True


def RemoveNull(matrix):
    """Function that removes rows consisting of just zeros"""

    null_rows = []

    for i in range(r := len(matrix)):
        null = 0
        for j in range(len(matrix[0])):
            if matrix[i][j] != 0:
                break
            null += 1

        if null != 0:
            null_rows.append(i)

    clean_matrix = []
    for i in range(r):
        if i not in null_rows:
            clean_matrix.append(matrix[i])

    return numpy.array(clean_matrix)
