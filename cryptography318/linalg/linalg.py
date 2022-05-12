from functools import reduce
from sympy import Symbol, im, solve
from numbers import Integral
from typing import Iterable
from cryptography318.utils.utils import where, shape
import numpy as np


def print_matrix(m, ret=False, rnd=3):
    """
    Prints matrix in readable format.

    :param m: Matrix instance
    :param ret: if function should return string instead of printing
    :param rnd: decimal places to round floats
    :return: Matrix as string if ``ret``, otherwise None
    """
    str_array = []
    max_len = 0
    shape = m.shape
    if len(shape) == 1:
        formatted = "[]"
    else:
        for i in range(shape[0]):
            str_array.append([])
            for j in range(shape[1]):
                v = m[i][j]

                # If not integral, make sure we round to avoid long floats
                if not isinstance(v, Integral):
                    v = round(v, rnd)
                s = str(v)
                max_len = max(len(s), max_len)
                str_array[i].append(s)

        # Make sure padding leaves room for at least 1 space AND is odd so that there is a center
        padding = (max_len + 1) | 1
        formatted = "["

        for i in range(shape[0]):
            # If not the first row, add a single space
            if i != 0:
                formatted += " "
            formatted += "["

            for j in range(len(str_array[0])):
                e = str_array[i][j]
                d = padding - len(e)
                pad_left = d // 2
                pad_right = d - pad_left
                formatted += pad_left * " " + e + " " * pad_right

            formatted += "]"

            # If it's not last row, add newline
            if i != shape[0] - 1:
                formatted += "\n"
        formatted += "]"

    if ret:
        return formatted
    else:
        print(formatted)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b, strict=True))


def transpose(a):
    return list(map(lambda i: list(map(lambda r: r[i], a)), range(len(a[0]))))


def matmul(a, b):
    T = transpose(b)
    return [[dot(x, y) for y in T] for x in a]


def flatten(a):
    if isinstance(a, Iterable):
        return reduce(lambda r1, r2: flatten(r1) + flatten(r2), a, [])
    else:
        return [a]


def matrix_slice(a, index):
    """
    Slices matrix by column index, returning the sliced matrix as a tuple
    with the first argument being the matrix from column 0 to index exclusive,
    and the second being matrix from index to last column.

    :param a: matrix
    :param index: column index of slice
    :return: matrix [0, index), matrix [index, )
    """
    left, right = [], []
    for row in a:
        left.append(row[:index])
        right.append(row[index:])
    return left, right


def matrix_copy(a):
    copy = []
    for row in a:
        copy.append(row[:])
    return copy


def matrix_equals(a, b):
    if len(a) == len(b):
        for x, y in zip(a, b):
            if len(x) != len(y) or any(i != j for i, j in zip(x, y)):
                return False
        return True
    else:
        return False


def make_pivot(a, index=None):
    if index is None:
        index = where(a)[0]
    return list(map(lambda n: n / a[index], a))


def row_reduce(a, row, col):
    w = len(a[row])
    for i in range(len(a)):
        if i != row:
            k = a[i][col]
            a[i] = [a[i][j] - a[row][j] * k for j in range(w)]


def identity_matrix(size):
    matrix = [[]] * size
    for i in range(size):
        matrix[i] = [0 if j != i else 1 for j in range(size)]
    return matrix


def ref(a, offset=0):
    pivot_row = 0  # first pivot belongs in first row

    w = len(a[0]) - offset
    h = len(a)

    array = matrix_copy(a)

    for j in range(w):

        # start at looking for pivot after previous pivot row
        for i in range(pivot_row, h):

            # if non-zero element, this row can become pivot row
            if array[i][j] != 0:

                if i > pivot_row:  # if pivot row not already in correct position, swap
                    array[i], array[pivot_row] = array[pivot_row], array[i]

                row_reduce(array, pivot_row, j)  # row reduce everything else
                pivot_row += 1
                break

    return array


def rref(a, offset=0):
    """
    Computes the reduced row-echelon form of the matrix. If offset is provided,
    computes the RREF ignoring the last offset columns.

    :param a: matrix
    :param offset: column index offset from last column
    :return: matrix in RREF
    """
    pivot_row = 0  # first pivot belongs in first row

    w = len(a[0]) - offset
    h = len(a)

    array = matrix_copy(a)

    for j in range(w):

        # start at looking for pivot after previous pivot row
        for i in range(pivot_row, h):

            # if non-zero element, this row can become pivot row
            if array[i][j] != 0:

                # make j'th element the pivot, reducing rest of row as well
                array[i] = make_pivot(array[i], j)
                if i > pivot_row:  # if pivot row not already in correct position, swap
                    array[i], array[pivot_row] = array[pivot_row][:], array[i][:]

                row_reduce(array, pivot_row, j)  # row reduce everything else
                pivot_row += 1
                break

    return array


def kernel(a):
    h = len(a)  # get number of rows
    w = len(a[0])
    array = list(list(r) for r in a)
    for j in range(w):  # this loop appends identity matrix to bottom of instance matrix
        row = [0] * w
        row[j] = 1
        array.append(row)

    array = rref(transpose(array), w)
    array, kern = matrix_slice(array, h)  # separates original matrix from now modified identity matrix
    basis = []

    for i, row in enumerate(array):
        if not any(row):  # all null rows in original matrix correspond to basis vector for kernel
            basis.append(kern[i])

    # Returns list of basis vectors for kernel i.e. rows are the vectors, not columns
    return basis


def binary_kernel(a):
    h = len(a)  # get number of rows
    w = len(a[0])
    array = matrix_copy(a)
    for j in range(w):  # this loop appends identity matrix to bottom of instance matrix
        row = [0] * w
        row[j] = 1
        array.append(row)

    array = transpose(array)

    pivot_row = 0  # first pivot belongs in first row

    for j in range(h):  # iterate only over current matrix, not attached identity matri

        # start at looking for pivot after previous pivot row
        for i in range(pivot_row, w):

            # if non-zero element, this row can become pivot row
            if array[i][j]:

                if i > pivot_row:  # if pivot row not already in correct position, swap
                    array[i], array[pivot_row] = array[pivot_row][:], array[i][:]

                for k in range(w):
                    if k != pivot_row and array[k][j]:
                        for m in range(h + w):
                            array[k][m] = (array[k][m] - array[pivot_row][m]) % 2
                pivot_row += 1

    array, kern = matrix_slice(array, h)  # separates original matrix from now modified identity matrix
    basis = []

    for i, row in enumerate(array):
        if not any(row):  # all null rows in original matrix correspond to basis vector for kernel
            basis.append(kern[i])

    # Returns list of basis vectors for kernel i.e. rows are the vectors, not columns
    return basis


def is_square(a):
    s = shape(a)
    return len(s) == 2 and s[0] == s[1]


def det(a):
    return float(np.linalg.det(np.array(a)))  # using numpy's det function for efficiency


def minor(a, index=None):
    """
    Computes the determinant of instance matrix if no index given. If index given,
    computes the minor A1,index (referring to the resulting matrix after removing
    row 1 and column index) multiplied against the value of A[0][index] with the
    correct sign (negative if index is even [when start counting at 1] otherwise
    positive). Use det(A) for calculating determinant for efficiency, unless specific
    minors or solving for sympy.Symbol is required.

    Examples
    --------
    >>> x = Symbol('x')
    >>> A = [[-1 - x, -3, 1], [3, 3 - x, 1], [3, 0, 4 - x]]
    >>> minor(A)
    -6*x + (3 - x)*(4 - x)*(-x - 1) + 18

    >>> A = [[-1, -3, 1], [3, 3, 1], [3, 0, 4]]
    >>> minor(A)
    6

    >>> A = [[-1, -3, 1], [3, 3, 1], [3, 0, 4]]
    >>> minor(A, 1)
    27
    >>> B = [[3, 1], [3, 4]]

    Notes
    -----
    in the first example, 'x' is a sympy.Symbol and the solution given is solvable using sympy.solve()

    in the last two examples, A.minor(1) computes the A[0, 1] * determinant of A1,2 (A with row 1 and column 2
    removed, start counting at 1) while B.det() computes the determinant of B, which is A1,2. Since the sign
    of the minor alternates, A.minor(1) returns -3 * -1 * det(A1,2) = -3 * -1 * B.det() = 27
    """

    if len(s := shape(a)) != 2 and len(set(s)) != 1:
        raise ValueError(f"{minor.__name__} incompatible with matrix of shape {s}")
    elif len(a) == 2:
        return a[0][0] * a[1][1] - a[1][0] * a[0][1]
    elif index is not None:
        return pow(-1, index % 2) * a[0][index] * minor([[e for i, e in enumerate(row) if i != index] for row in a[1:]])
    else:
        return sum(minor(a, j) for j in range(len(a[0])))


def char_poly(a, sym=None):
    """Computes the characteristic polynomial of a square matrix. Analogous
    to A.minor() for A = instance-matrix - Identity * x, for some variable x
    [typically x = sympy.Symbol('x')]."""

    if len(s := shape(a)) != 2 and len(set(s)) != 1:
        raise ValueError("Matrix must be square!")
    elif sym is None:
        sym = Symbol("x")
    elif not isinstance(sym, Symbol):
        sym = Symbol(sym)

    h = len(a)
    return minor([[a[i][j] - sym if i == j else a[i][j] for j in range(h)] for i in range(h)])


def eigvals(a):
    """Finds all eigenvalues for square matrix. Solves equation det(A - xI) = 0 with A
    equal to the instance matrix, I the identity, for x. The set of all solutions for x
    is analogous to the set of eigenvalues for A.

    :return: list of real or imaginary eigenvalues

    Examples
    --------
    >>> A = [[-1, -3, 1], [3, 3, 1], [3, 0, 4]]
    >>> eigvals(A)
    [1.0, 2.0, 3.0]

    >>> A = [[-1, -3, 1], [3, 3, 1], [3, 0, 3]]
    >>> eigvals(A)
    [0, (2.5-1.6583123951777j), (2.5+1.6583123951777j)]
    """

    if len(s := shape(a)) != 2 and len(set(s)) != 1:
        raise AttributeError("Matrix must be square")

    x = Symbol('x')
    return [float(e) if im(e) == 0 else complex(e) for e in solve(char_poly(a, x), x)]


def eigvec(a, values=None):
    if len(s := shape(a)) != 2 and len(set(s)) != 1:
        raise AttributeError("Matrix must be square")
    elif values is None:
        values = eigvals(a)

    vectors = {}
    for v in values:
        mat = [[e - v if i == j else e for i, e in enumerate(row)] for j, row in enumerate(a)]
        vec = kernel(mat)
        if len(vec) == 1:
            vectors[v] = vec[0]
        else:
            vectors[v] = vec
    return vectors
