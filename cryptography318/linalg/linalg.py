import numpy as np


def kernel_gf2(a):
    """
    Computes the null space of a matrix in GF(2)
    :param a: matrix
    :return: basis of null space as rows
    """

    # dtype=object, so we don't lose any precision if being passed large Python integers
    a = np.array(a, dtype=object) & 1
    h, w = a.shape

    # Append identity matrix to bottom
    array = np.append(a, np.eye(w), axis=0).transpose().astype(np.int8)

    pivot_row = 0  # first pivot belongs in first row

    for j in range(h):  # iterate only over current matrix, not attached identity matrix

        # start at looking for pivot after previous pivot row
        for i in range(pivot_row, w):

            # if non-zero element, this row can become pivot row
            if array[i, j] == 1:

                if i > pivot_row:  # if pivot row not already in correct position, swap
                    array[[i, pivot_row]] = array[[pivot_row, i]]

                f = np.eye(w, dtype=np.int8)
                f[:, pivot_row] = array[:, j]
                array = (f @ array) & 1
                pivot_row += 1
                break

    # What we are left with in 'array' is a matrix w by h in rref on the left and
    # a matrix w by w on the right containing our potential kernel.

    left, _, right = np.split(array, [h, h], axis=1)

    # Our kernel is all the rows on the right where the row on the left is all 0's
    # and since our matrix on the left is in rref, as soon as we find one null row
    # all the rows after that cannot contain pivots either
    i = 0
    while i < w and any(left[i]):
        i += 1

    # Now that we have found `i` being the first index of a null row, return the rest
    return right[i:]


def kernel(a):
    """
    Computes null space of matrix
    :param a: matrix
    :return: basis of null space as rows
    """
    array = np.array(a)
    if len(array.shape) != 2:
        raise ValueError(f"{kernel.__name__} can only be performed on 2-dimensional arrays")

    u, s, vh = np.linalg.svd(array)
    rank, = s.shape
    return vh[rank + 1:]


def rref(a, offset=0, rtol=1e-05, atol=1e-08):
    """
    Computes the reduced row-echelon form of the matrix. If offset is provided,
    computes the RREF ignoring the last offset columns.

    :param a: matrix
    :param offset: column index offset from last column
    :param rtol:
    :param atol:
    :return: matrix in RREF
    """

    shape = np.shape(a)
    if len(shape) != 2:
        raise ValueError(f"{rref.__name__} can only be performed on 2-dimensional arrays")

    h, w = shape
    pivot_row = 0  # first pivot belongs in first row

    array = np.array(a)

    for j in range(w - offset):

        # start at looking for pivot after previous pivot row
        for i in range(pivot_row, h):

            # if non-zero element, this row can become pivot row
            if not np.isclose(array[i, j], 0, rtol=rtol, atol=atol):

                # make j'th element the pivot, reducing rest of row as well
                array[i] /= array[i, j]
                if i > pivot_row:  # if pivot row not already in correct position, swap
                    array[[i, pivot_row]] = array[[pivot_row, i]]

                # This process row reduces
                f = np.eye(h)
                b = -array[:, j]  # Set the active column to negative values for each row
                b[pivot_row] = 1  # And positive 1 for the pivot row
                f[:, pivot_row] = b
                array = f @ array  # Matrix multiplication means each row is subtracted by the pivot row

                pivot_row += 1
                break

    return array


"""
Functions below here are not declared in __init__ as they are untested and not optimized.
TODO: Incorporate these
"""


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
    >>> import sympy
    >>> x = sympy.Symbol('x')
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

    if len(s := np.shape(a)) != 2 and len(set(s)) != 1:
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

    from sympy import Symbol

    if len(s := np.shape(a)) != 2 and len(set(s)) != 1:
        raise ValueError("Matrix must be square!")
    elif sym is None:
        sym = Symbol("x")
    elif not isinstance(sym, Symbol):
        sym = Symbol(sym)

    h = len(a)
    return minor([[a[i][j] - sym if i == j else a[i][j] for j in range(h)] for i in range(h)])


def eigvals(a):
    """Finds all eigenvalues for square matrix. Solves equation det(A - xI) = 0 with A
    equal to the instance matrix, `I` the identity, for x. The set of all solutions for x
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

    from sympy import Symbol, solve, im

    if len(s := np.shape(a)) != 2 and len(set(s)) != 1:
        raise AttributeError("Matrix must be square")

    x = Symbol('x')
    return [float(e) if im(e) == 0 else complex(e) for e in solve(char_poly(a, x), x)]


def eigvec(a, values=None):
    if len(s := np.shape(a)) != 2 and len(set(s)) != 1:
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
