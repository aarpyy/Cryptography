import sys

from typing import Callable
from functools import reduce
from typing import Sequence, Iterable
from numpy import array as np_array
from numpy.linalg import det as np_det
from numpy import where
from sympy import Symbol, im, solve


def dot_3_10(a: Sequence, b: Sequence):
    return sum(x * y for x, y in zip(a, b, strict=True))


def dot_3_9(a: Sequence, b: Sequence):
    if len(a) == len(b):
        return sum(x * y for x, y in zip(a, b))
    else:
        raise ValueError(f"Lengths {len(a)} and {len(b)} differ")


if sys.version.startswith("3.10"):
    dot = dot_3_10      # type: Callable
else:
    dot = dot_3_9       # type: Callable


def transpose(a: Sequence[Sequence]):
    return list(map(lambda i: list(map(lambda r: r[i], a)), range(len(a[0]))))


def vec_matmul(a: Sequence, b: Sequence[Sequence]):
    return [dot(a, x) for x in transpose(b)]


def matmul(a: Sequence[Sequence], b: Sequence[Sequence]):
    return [vec_matmul(x, b) for x in a]


def flatten(a: Iterable):
    return [*reduce(lambda r1, r2: list(r1) + list(r2), a)]


def matrix_slice(a, index=0):
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


def make_pivot(a, index=None):
    if index is None:
        index = where(a)[0][0]
    return list(map(lambda n: n / a[index], a))


def row_reduce(a, row, col):
    h = len(a)
    w = len(a[row])
    for i in range(h):
        if i != row:
            for j in range(w):
                a[i][j] -= a[row][j] * a[i][col]


def identity_matrix(size):
    matrix = [[]] * size
    for i in range(size):
        matrix[i] = [0 if j != i else 1 for j in range(size)]
    return matrix


def rref(a):

    pivot_row = 0  # first pivot belongs in first row

    w = len(a[0])
    h = len(a)

    array = matrix_copy(a)

    for j in range(w):

        # start at looking for pivot after previous pivot row
        for i in range(pivot_row, h):

            # if non-zero element, this row can become pivot row
            if array[i][j] != 0:

                # make j'th element the pivot, reducing rest of row as well
                make_pivot(array[i], j)
                if i > pivot_row:  # if pivot row not already in correct position, swap
                    array[i], array[pivot_row] = array[pivot_row][:], array[i][:]

                row_reduce(array, pivot_row, j)  # row reduce everything else
                pivot_row += 1
                break

    return array


def kernel(a):
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
        for i in range(pivot_row, len(array)):

            # if non-zero element, this row can become pivot row
            if array[i][j] != 0:

                # make j'th element the pivot, reducing rest of row as well
                make_pivot(array[i], j)
                if i > pivot_row:  # if pivot row not already in correct position, swap
                    array[i], array[pivot_row] = array[pivot_row][:], array[i][:]

                row_reduce(array, pivot_row, j)  # row reduce everything else
                pivot_row += 1

    array, kern = matrix_slice(array, h)  # separates original matrix from now modified identity matrix
    basis = []

    for i, row in enumerate(array):
        if not any(row):    # all null rows in original matrix correspond to basis vector for kernel
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
        if not any(row):    # all null rows in original matrix correspond to basis vector for kernel
            basis.append(kern[i])

    # Returns list of basis vectors for kernel i.e. rows are the vectors, not columns
    return basis


def det(a):
    return float(np_det(np_array(a)))  # using numpy's det function for efficiency


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

    h = len(a)
    if h == 2:
        return a[0][0] * a[1][1] - a[1][0] * a[0][1]

    w = len(a[0])
    if index is not None:
        m = -a[0][index] if index % 2 else a[0][index]
        a = a[1:]
        for i in range(h):
            a[i] = [m * a[i][j] for j in range(w) if j != index]
        return minor(a)

    det = 0
    sign = 1
    for j in range(w):

        # Increment det with the product of the first row, j'th column element and the
        # minor of a with the first row and j'th column missing
        b = a[1:]
        for i in range(h):
            b[i] = [a[0][j] * b[i][k] for k in range(w) if k != j]
        det += sign * minor(b)
        sign *= -1

    return det


def char_poly(a, sym='x'):
    """Computes the characteristic polynomial of a square matrix. Analogous
    to A.minor() for A = instance-matrix - Identity * x, for some variable x
    [typically x = sympy.Symbol('x')]."""

    if len(a) != len(a[0]):
        raise ValueError("Matrix must be square!")

    if not isinstance(sym, Symbol):
        sym = Symbol(sym)

    h = len(a)
    return [[a[i][j] - sym if i == j else a[i][j] for j in range(h)] for i in range(h)]


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

    if len(a) != len(a[0]):
        raise AttributeError("Matrix must be square")

    x = Symbol('x')
    return [float(e) if im(e) == 0 else complex(e) for e in solve(char_poly(a, x), x)]


def eigvec(a, values=None):
    if len(a) != len(a[0]):
        raise AttributeError("Matrix must be square")

    if values is None:
        values = eigvals(a)

    vectors = []
    for v in values:
        mat = [[e - v for e in row] for row in a]
        kern = kernel(mat)
        for vec in kern:
            vectors.append(vec)
    return vectors

