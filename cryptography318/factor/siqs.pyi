from typing import Optional, Iterable, Sequence
from random import Random

MIN_A_FACTOR: int
MAX_A_FACTOR: int
MIN_N_FACTORS: int
TRIALS_A: int
TRIAL_ERROR_MARGIN: int
REQUIRED_RELATIONS_RATIO: float
TRIALS_LINALG: int
relations_found: int
min_sieve: int
loud_print: bool
a: int
b: int
B: list[int]
primes: list[int]
t_sqrt: list[int]
log_p: list[int]
factor_base: list[int]
soln1: list[int]
soln2: list[int]
B_ainv_2: list[list[int]]
a_factors: set[int]
a_non_factors: set[int]
smooth_u: list[list[int]]
smooth_t: list[int]
rand: Random

class QSPoly:
    args: list[int]
    def __new__(cls, *args: int) -> QSPoly: ...
    def __call__(self, x: int) -> int: ...

def l_print(*args, **kwargs) -> None: ...
def choose_f(digits: int) -> int: ...
def choose_m(digits: int) -> int: ...
def init_siqs(n, *, fp: Optional[str] = ...) -> None: ...
def smooth_a(n: int, m: int) -> None: ...
def first_poly(n: int, m: int) -> tuple[QSPoly, QSPoly]: ...
def next_poly(i: int, n: int) -> tuple[QSPoly, QSPoly]: ...
def sieve(m: int) -> list[int]: ...
def trial_division(sieve_array: Iterable[int], m: int, g: QSPoly, h: QSPoly) -> None: ...
def vec_matmul_T(vector: Sequence[int], matrix: Sequence[Sequence[int]]) -> list[int]: ...
def solve_matrix(n: int) -> Optional[int]: ...
def siqs(n: int, *, fp: str = ..., loud: bool = ...) -> Optional[int]: ...
