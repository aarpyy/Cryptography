from .linalg import (
    rref, kernel, binary_kernel, minor, det, eigvec, eigvals, char_poly, flatten, matmul, transpose, dot, matrix_equals
)
from .linalg_np import kernel_gf2 as binary_kernel_np

__all__ = [
    "rref", "kernel", "binary_kernel", "minor", "det", "eigvals", "eigvec", "char_poly", "flatten",
    "matmul", "transpose", "dot", "matrix_equals", "binary_kernel_np"
]
