from typing import Optional, Iterable, Sequence

def smooth_factor(n: int, factors: Sequence[int]) -> Optional[list[int]]: ...
def eval_power(exp: Iterable[int], primes: Iterable[int]) -> int: ...
def from_base(lst: Iterable[int], base: int) -> int: ...
def extended_gcd(*args: int) -> tuple[int]: ...
def n_digits(n: int) -> int: ...