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


def ArrayToList(matrix):
    lst = []
    for i in range(len(matrix)):
        lst.append([])
        for j in range(len(matrix[0])):
            lst[i].append(matrix[i][j])
    return lst


def RandomMatrix(rows=None, cols=None):
    if rows is None:
        rows = random.randrange(1, 10)
    if cols is None:
        cols = random.randrange(1, 10)

    return MakeMatrix(rows, cols, rand=True)


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


def MatrixFloat(matrix, bits=64):
    float_types = {16: numpy.float16, 32: numpy.float32, 64: numpy.float64}
    if bits in float_types:
        return matrix.astype(float_types[bits])
    return matrix.astype(float)


def ResetType(matrix):
    """Takes as input a matrix of unknown number type and returns a matrix
    consisting of type integers if possible, type numpy.float16 if not."""

    matrix = ArrayToList(matrix.astype(numpy.float16))
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            e = str(matrix[i][j]).split(".")
            if len(e) > 1 and e[1] == '0':
                matrix[i][j] = int(matrix[i][j])
    return numpy.array(matrix)


def augmentMatrix(matrix, solutions):
    """Takes as input a coefficient matrix and a column matrix of solutions and
    returns one augmented coefficient matrix."""

    m = ArrayToList(matrix)
    for i in range(len(m)):
        m[i].append(solutions[i][0])
    return numpy.array(m)


def separateMatrix(matrix):
    """Takes as input an augmented coefficient matrix and returns two separate instances
    of the original coefficient matrix with a separate column matrix as the constants."""

    matrixM, sol = [], []
    c = len(matrix[0]) - 1
    for i in range(len(matrix)):
        matrixM.append(list(matrix[i][:c]))
        sol.append([matrix[i][c]][:])
    return numpy.array(matrixM), numpy.array(sol)


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


def Transpose(matrix):
    matrixT = []
    for j in range(len(matrix[0])):
        matrixT.append([])
        for i in range(len(matrix)):
            matrixT[j].append(matrix[i][j])
    return numpy.array(matrixT)
