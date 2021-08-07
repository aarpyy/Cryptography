import numpy
from .matrix_deprecated import augmentMatrix, MatrixFloat, ResetType, MakeMatrix, ArrayToList, MultiplyMatrix


def CopyMatrix(matrix):
    new_mat = []
    for i in range(len(matrix)):
        new_mat.append([])
        for j in range(len(matrix[0])):
            new_mat[i].append(matrix[i][j])

    return numpy.array(new_mat)


def RemoveNull(matrix):
    """Function that removes rows consisting of just zeros"""

    matrix = CopyMatrix(matrix)
    null_rows = []
    for i in range(c := len(matrix)):
        null = 0
        for j in range(r := len(matrix[0])):
            if matrix[i][j] != 0:
                break
            null += 1

        if null == r:
            null_rows.append(i)

    clean_matrix = []
    for i in range(c):
        if i not in null_rows:
            clean_matrix.append(matrix[i][:])

    return numpy.array(clean_matrix)


def ChangeBasis(matrix, basis):
    """Function that takes as input a matrix in a standard basis and returns
    the matrix in the given basis"""

    return MultiplyMatrix(InvertMatrix(basis), matrix)


def ChangeBasisMap(lin_map, in_basis=None, out_basis=None):
    """Function that takes as input a linear map and at least one basis and outputs
    the linear map with a change of basis.

    Linear map: an m x n matrix mapping from the standard basis vectors in m-dimensions to
    the standard basis vectors in n-dimensions.
    Notation: [S]E->E, where S, E are the linear map and standard basis vectors respectively.

    In-Basis: an m x m matrix that forms a basis in m-dimensions. This matrix, when applied
    to the linear map, gives the linear map such that the map takes as input matrices in
    this basis, and outputs them in the standard basis of n-dimensions.
    Notation: [S]B->E, where B is the In-Basis.

    Out-Basis: an n x n matrix that forms a basis in n-dimensions. This matrix, when applied
    to the linear map, gives the linear map such that the map takes as input matrices in
    the standard basis of m-dimensions and outputs a matrix in the form of the new basis.
    Notation: [S]E->C, where C is the Out-Basis"""

    if in_basis is None and out_basis is None:
        raise ValueError("Make sure to enter at least one basis")

    inputs = [lin_map, in_basis, out_basis]
    for e in inputs:
        if e is not None and type(e) is not numpy.ndarray:
            raise ValueError("All input must be type numpy.ndarray")

    if in_basis is None:
        if len(out_basis) != len(out_basis[0]):
            raise ValueError("Basis must be square matrix")
        return MultiplyMatrix(InvertMatrix(out_basis), lin_map)

    if out_basis is None:
        if len(in_basis) != len(in_basis[0]):
            raise ValueError("Basis must be square matrix")
        return MultiplyMatrix(lin_map, in_basis)

    if len(in_basis) != len(in_basis[0]) or len(out_basis) != len(out_basis[0]):
        raise ValueError("Bases must be square matrix")
    return MultiplyMatrix(InvertMatrix(out_basis), MultiplyMatrix(lin_map, in_basis))


def InvertMatrix(matrix=None):
    if matrix is None:
        size = int(input("Enter matrix dimensions nxn: "))
        matrix = MakeMatrix(size, size)
    if len(matrix) != len(matrix[0]):
        raise ValueError("InvertMatrix only works with square matrices")
    return ResetType(numpy.linalg.inv(CopyMatrix(matrix)))


def _SortRREF_old(matrix):
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


def _RREF_old(matrix=None):
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

    return ResetType(_SortRREF_old(matrixM))


def RREF(matrix=None):
    if matrix is None:
        print("Enter dimensions for Matrix: ")
        rows = int(input("Rows: "))
        cols = int(input("Cols: "))
        matrix = MakeMatrix(rows, cols)
    m = ArrayToList(MatrixFloat(CopyMatrix(matrix)))

    pivot_row = -1
    # iterates across columns
    for j in range(row_len := len(m[0])):
        # iterates down columns
        pivot = False
        for i in range(col_len := len(m)):
            e = m[i][j]
            # if found a non-zero entry that could be a pivot, and is not already a row with a pivot, make it a pivot
            if e != 0 and i > pivot_row:
                # shifts entire row by factor that makes first entry = 1, therefore is pivot
                # print(f"e = {e}, m[i] before shift: {m[i]}")
                # print(shift)
                for k in range(row_len):
                    m[i][k] /= e
                # print(m[i])
                # pivot row increases from where it last was, so that new row will be one below last pivot row
                pivot_row += 1
                # if pivot row isn't current one, swap so that row with pivot comes directly after previous pivot row
                if pivot_row != i:
                    m[i][:], m[pivot_row][:] = m[pivot_row][:], m[i][:]
                # this row is a pivot
                pivot = True
                break
        if pivot:
            # iterate down through matrix, removing all non-zero entries from column with pivot
            for k in range(col_len):
                e = m[k][j]
                if k != pivot_row and e != 0:
                    for l in range(row_len):
                        # here, e represents the number of pivot row's needed to be removed to make i'th row have
                        # a zero entry in this column, ex. pivot row has 1 in column, i'th row as 3, removing 3 of
                        # pivot row will make i'th row have 0 in column
                        m[k][l] -= m[pivot_row][l] * e

    return ResetType(numpy.array(m))


def REF(matrix=None):
    if matrix is None:
        print("Enter dimensions for Matrix: ")
        rows = int(input("Rows: "))
        cols = int(input("Cols: "))
        matrix = MakeMatrix(rows, cols)
    m = ArrayToList(MatrixFloat(CopyMatrix(matrix)))

    pivot_row = -1
    # iterates across columns
    for j in range(row_len := len(m[0])):
        # iterates down columns
        pivot = False
        for i in range(col_len := len(m)):
            e = m[i][j]
            # if found a non-zero entry that could be a pivot, and is not already a row with a pivot, make it a pivot
            if e != 0 and i > pivot_row:
                # shifts entire row by factor that makes first entry = 1, therefore is pivot
                # print(f"e = {e}, m[i] before shift: {m[i]}")
                # print(shift)
                for k in range(row_len):
                    m[i][k] /= e
                # print(m[i])
                # pivot row increases from where it last was, so that new row will be one below last pivot row
                pivot_row += 1
                # if pivot row isn't current one, swap so that row with pivot comes directly after previous pivot row
                if pivot_row != i:
                    m[i][:], m[pivot_row][:] = m[pivot_row][:], m[i][:]
                # this row is a pivot
                pivot = True
                break
        if pivot:
            # iterate down through matrix, removing all non-zero entries from column with pivot
            for k in range(col_len):
                e = m[k][j]
                if k > pivot_row and e != 0:
                    for l in range(row_len):
                        # here, e represents the number of pivot row's needed to be removed to make i'th row have
                        # a zero entry in this column, ex. pivot row has 1 in column, i'th row as 3, removing 3 of
                        # pivot row will make i'th row have 0 in column
                        m[k][l] -= m[pivot_row][l] * e

    return ResetType(numpy.array(m))


def augRREF(matrix=None):
    if matrix is None:
        print("Enter dimensions for Matrix: ")
        rows = int(input("Rows: "))
        cols = int(input("Cols: "))
        matrix = MakeMatrix(rows, cols)
    m = ArrayToList(MatrixFloat(CopyMatrix(matrix)))

    pivot_row = -1
    # iterates across all columns except last
    for j in range((row_len := len(m[0])) - 1):
        # iterates down columns
        pivot = False
        for i in range(col_len := len(m)):
            e = m[i][j]
            # if found a non-zero entry that could be a pivot, and is not already a row with a pivot, make it a pivot
            if e != 0 and i > pivot_row:
                # shifts entire row by factor that makes first entry = 1, therefore is pivot
                # print(f"e = {e}, m[i] before shift: {m[i]}")
                # print(shift)
                for k in range(row_len):
                    m[i][k] /= e
                # print(m[i])
                # pivot row increases from where it last was, so that new row will be one below last pivot row
                pivot_row += 1
                # if pivot row isn't current one, swap so that row with pivot comes directly after previous pivot row
                if pivot_row != i:
                    m[i][:], m[pivot_row][:] = m[pivot_row][:], m[i][:]
                # this row is a pivot
                pivot = True
                break
        if pivot:
            # iterate down through matrix, removing all non-zero entries from column with pivot
            for k in range(col_len):
                e = m[k][j]
                if k != pivot_row and e != 0:
                    for l in range(row_len):
                        # here, e represents the number of pivot row's needed to be removed to make i'th row have
                        # a zero entry in this column, ex. pivot row has 1 in column, i'th row as 3, removing 3 of
                        # pivot row will make i'th row have 0 in column
                        m[k][l] -= m[pivot_row][l] * e

    return ResetType(numpy.array(m))


def _augRREF_old(matrix, solutions=None):
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

    return ResetType(_SortRREF_old(matrixM))


def IsRREF1(matrix):
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
        if pivots[c] > pivots[c + 1]:
            return False
    return True


def IsRREF(matrix):
    """Function that returns True if a matrix is in reduced-row echelon form, and False otherwise"""

    pivot_col, pivot_row = [], []
    for j in range(len(matrix[0])):
        for i in range(len(matrix) - 1, -1, -1):
            e = matrix[i][j]
            if e != 0 and j in pivot_col:
                return False
            if e == 1 and i not in pivot_row:
                pivot_col.append(j)
                pivot_row.append(i)
    return True


def augIsRREF(matrix):
    """Function that returns True if a matrix is in reduced-row echelon form, and False otherwise"""

    pivot_col, pivot_row = [], []
    for j in range(len(matrix[0]) - 1):
        for i in range(len(matrix) - 1, -1, -1):
            e = matrix[i][j]
            if e != 0 and j in pivot_col:
                return False
            if e == 1 and i not in pivot_row:
                pivot_col.append(j)
                pivot_row.append(i)
    return True


def IsConsistent(matrix):
    """Function that takes as input an augmented coefficient matrix and returns
        True if the matrix is consistent, False if otherwise"""

    if not augIsRREF(matrix):
        matrix = augRREF(matrix)

    for i in range(len(matrix) - 1, -1, -1):
        for j in range(r := len(matrix[0])):
            e = matrix[i][j]
            if e == 1 and j < r - 1:
                return True
            if e != 0 and j == r - 1:
                return False
    return True


def IsSolvable(matrix, aug=False):
    if aug and not augIsRREF(matrix):
        matrix = augRREF(matrix)
    elif not aug and not IsRREF(matrix):
        matrix = RREF(matrix)
    m = RemoveNull(matrix)
    if aug:
        if len(m) == len(m[0]) - 1:
            return True
        return False
    if len(m) == len(m[0]):
        return True
    return False


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
            print("this false, 422")
            return False

        for x in _vars:
            _vars[x] = matrix[0][r]
        return _vars

    result = _solveSystem(matrix[1:], sol=None, aug=True)
    if result is False:
        print("this false, 431")
        return False

    unknowns = []
    for i in range(r):
        if matrix[0][i] != 0:
            if i not in result:
                unknowns.append(i)

    if len(unknowns) > 1:
        print("this false, 441")
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


def IsSolvable1(matrix, aug=False):
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
            if abs(round(mat1[i][j]) - round(mat2[i][j])) > 1:
                return False
    return True
