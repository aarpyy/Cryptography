import numpy as np


def kernel_gf2(a: np.ndarray):
    """
    See /Users/andrew/Downloads/1006.1744.pdf for potentially more efficient replacements
    :param a:
    :return:
    """
    h, w = a.shape

    # Append identity matrix to bottom
    array = np.append(a, np.eye(w, dtype=np.int8), axis=0).transpose()

    pivot_row = 0  # first pivot belongs in first row

    for j in range(h):  # iterate only over current matrix, not attached identity matrix

        # start at looking for pivot after previous pivot row
        for i in range(pivot_row, w):

            # if non-zero element, this row can become pivot row
            if array[i][j] == 1:

                if i > pivot_row:  # if pivot row not already in correct position, swap
                    array[[i, pivot_row]] = array[[pivot_row, i]]

                for k in range(w):
                    if k != pivot_row and array[k][j]:
                        array[k][:] = (array[k] - array[pivot_row]) % 2
                pivot_row += 1
                break

    i = 0
    # What we are left with in 'array' is a matrix w by h in rref on the left and
    # a matrix w by w on the right containing our potential kernel.

    # Our kernel is all of the rows on the right where the row on the left is all 0's
    # and since our matrix on the left is in rref, as soon as we find one null row
    # all of the rows after that cannot contain pivots either
    while any(array[i, :h]):
        i += 1

    # Now that we have found i being the first index of a null row, return the rest
    return array[i:, h:]
