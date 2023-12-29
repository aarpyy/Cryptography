from collections.abc import Callable
from os import PathLike
from random import Random


class FactorBaseItem:
    prime: int
    soln1: int | None
    soln2: int | None
    log_p: int
    B_ainv_2: list[int] | None
    t_sqrt: int
    a_inv: int | None

    def __init__(self, prime, soln1, soln2, log_p, B_ainv_2, t_sqrt, a_inv) -> None: ...


MIN_A_FACTOR: int
MAX_A_FACTOR: int
MIN_N_FACTORS: int
TRIALS_A: int
TRIAL_ERROR_MARGIN: int
REQUIRED_RELATIONS_RATIO: float
TRIALS_LINALG: int
rand: Random


class QSPoly(Callable[[int], int]):
    args: tuple[int, ...]
    a: int
    b: int

    def __new__(cls, *args, a: int | None = ..., b: int | None = ...): ...

    def __call__(self, x: int) -> int: ...


def vec_matmul_T(vector, matrix): ...


class SIQS:
    verbose: bool
    n: int
    factor_base: list[FactorBaseItem]
    required_relations: int
    m: int
    smooth_u: list[list[int]]
    smooth_t: list[int]

    def __init__(self, n: int, F: int, *, fp: int | str | bytes | PathLike[str] | PathLike[bytes] | None = ...,
                 verbose: bool = ...) -> None: ...

    def print(self, *args, **kwargs) -> None: ...

    def smooth_a(self) -> tuple[int, set[int]]: ...

    def first_poly(self): ...

    def next_poly(self, i: int, g: QSPoly, B: list[int]): ...

    def sieve(self) -> list[int]: ...

    def smooth_factor(self, x: int) -> list[int] | None: ...

    def trial_division(self, sieve_array: list[int], g: QSPoly) -> int: ...

    def solve_matrix(self) -> int | None: ...

    def clear(self) -> None: ...


def qs(n: int, F: int = ..., *, fp: int | str | bytes | PathLike[str] | PathLike[bytes] | None = ..., verbose: bool = ...): ...
