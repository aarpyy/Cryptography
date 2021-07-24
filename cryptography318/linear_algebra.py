import numpy
from .matrix import augmentMatrix, MatrixFloat, ResetType, MakeMatrix, ArrayToList


def CopyMatrix(matrix):
    new_mat = []
    for i in range(len(matrix)):
        new_mat.append([])
        for j in range(len(matrix[0])):
            new_mat[i].append(matrix[i][j])

    return numpy.array(new_mat)


def InvertMatrix(matrix=None):
    if matrix is None:
        size = int(input("Enter matrix dimensions nxn: "))
        matrix = MakeMatrix(size, size)
    return numpy.linalg.inv(CopyMatrix(matrix))


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
    matrixM = MatrixFloat(CopyMatrix(matrix), 64)

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


def IsConsistent(matrix):
    """Function that takes as input an augmented coefficient matrix and returns
    True if the matrix is consistent, False if otherwise"""

    c = len(matrix[0]) - 1
    for i in range(len(matrix)):
        non_zero = 0
        for j in range(c):
            if matrix[i][j] != 0:
                non_zero += 1
        if non_zero == 0 and matrix[i][c] != 0:
            return False
    return True


def _solveSystem(matrix, sol=None, aug=False):
    """
    _solveSystem is a private helper function for public Solve that recursively
    determines solution for each variable in system of linear equations. Function
    calls itself recursively with input matrix being the given matrix minus the
    first row, until the input matrix is one row. Solving this row is then attempted,
    returning False if variable has no known solution, and a dictionary with the variable
    as the key and solution as the value otherwise. At each recursive step, the next
    variable and its solution are again added as the key and value in the dictionary.
    """
    if sol is None and not aug:
        sol = numpy.array([[0]] * len(matrix[0]))

    if not aug:
        matrix = augmentMatrix(matrix, sol)

    if not augIsRREF(matrix):
        matrix = augRREF(matrix)

    r = len(matrix[0]) - 1

    _vars = {}
    if len(matrix) == 1:
        for i in range(r):
            if matrix[0][i] != 0:
                _vars[i] = None
        if len(_vars) > 1:
            return False

        for x in _vars:
            _vars[x] = matrix[0][r]
        return _vars

    result = _solveSystem(matrix[1:], sol=None, aug=True)
    if result is False:
        return False

    unknowns = []
    for i in range(r):
        if matrix[0][i] != 0:
            if i not in result:
                unknowns.append(i)

    if len(unknowns) > 1:
        return False
    if len(unknowns) == 0:
        return result

    solution = 0
    index = None
    for i in unknowns:
        if i in result:
            solution += result[i] * matrix[0][1]
        else:
            index = i

    result[index] = matrix[0][r] - solution
    return result


def Solve(matrix, sol=None):
    """
    Solve is a function that attempts to solve a given system of linear equations

    Input is expected to be in the form:
        matrix->numpy.array, either an augmented coefficient matrix with sol=None or standard matrix
        sol->numpy.array, represents the set of solutions, if sol is None matrix is expected to be in
        the form of an augmented coefficient matrix
    """
    aug = True if sol is None else False
    if not IsSolvable(matrix, aug):
        return False

    # result could either be False, empty list, or list of solutions
    result = _solveSystem(MatrixFloat(matrix), sol, aug)
    if not result:
        return result

    solution = []
    for var in sorted(result):
        solution.append([result[var]])

    return ResetType(numpy.array(solution))


def IsSolvable(matrix, aug=False):
    """Determines if a given system of linear equations is solvable. First, matrix
     is checked to see if it is already in reduced-row echelon form, putting it in RREF if not.
     If matrix is an augmented coefficient (aug=True) then matrix is also checked for consistency.
     Then matrix is checked for the number of free variables, returning False if any are found."""

    if aug:
        if not augIsRREF(matrix):
            matrix = augRREF(matrix)
        if not IsConsistent(matrix):
            return False
    else:
        if not IsRREF(matrix):
            matrix = RREF(matrix)

    r = len(matrix[0])
    if aug:
        r -= 1

    for j in range(r):
        for i in range(len(matrix)):
            if matrix[i][j] not in [0, 1]:
                return False
    return True


def MatrixEquals(mat1, mat2):
    """Returns True if two given matrices are equal to a zero decimal point precision.
    All numbers are converted into integers for comparison, since any two matrices
    that are dis-similar by less than a whole integer value, they are considered equal."""

    if len(mat1) != len(mat2) or len(mat1[0]) != len(mat2[0]):
        return False

    for i in range(len(mat1)):
        for j in range(len(mat1[0])):
            if int(mat1[i][j]) != int(mat2[i][j]):
                return False
    return True
