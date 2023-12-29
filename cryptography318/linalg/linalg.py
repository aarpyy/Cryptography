import numpy as np


def kernel_gf2(a, dtype=object):
    """
    Computes the null space of a matrix in GF(2)
    :param a: matrix
    :param dtype: numpy dtype of matrix
    :return: basis of null space as rows
    """

    # dtype=object, so we don't lose any precision if being passed large Python integers
    a = np.array(a, dtype=dtype) & 1
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
