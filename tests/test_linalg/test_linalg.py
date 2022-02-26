from linalg import *


def test_matmul():
    a = [[1, 2, 3], [4, 5, 6]]
    b = [[1, 2], [1, 2], [1, 2]]
    assert matmul(a, b) == [[6, 12], [15, 30]]


def test_matrix_slice():
    a = [[1, 2, 3], [4, 5, 6]]
    left, right = matrix_slice(a, 1)
    assert left == [[1], [4]]
    assert right == [[2, 3], [5, 6]]


def test_binary_kernel():
    matrix = [[1, 0, 1, 0, 0, 0],
              [0, 1, 1, 0, 1, 0],
              [0, 0, 0, 1, 1, 1],
              [0, 0, 0, 0, 0, 0]]

    kern = binary_kernel(matrix)
    print(kern)
    for row in kern:
        assert all(dot(row, e) % 2 == 0 for e in matrix)


if __name__ == "__main__":
    test_matmul()

