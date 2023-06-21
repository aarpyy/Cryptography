from typing import Any, Iterable, Optional, Sequence


class Details:

    _details: dict[str, Any]
    def __init__(self, **kwargs) -> None: ...
    def __str__(self) -> str: ...
    def add_details(self, name: str, value: Any) -> None: ...
    def clear_details(self) -> None: ...

def binary_search(a: Sequence, key: Any, *, start: int = ..., end: int | None = ...,
                  exist: bool = ...) -> int | None: ...


def smooth_factor(n: int, factors: Sequence[int]) -> Optional[list[int]]: ...


def eval_power(exp: Iterable[int], primes: Iterable[int]) -> int: ...


def from_base(lst: Iterable[int], base: int) -> int: ...


def extended_gcd(*args: int) -> tuple[int]: ...


def n_digits(n: int) -> int: ...
