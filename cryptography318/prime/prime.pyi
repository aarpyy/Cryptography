from typing import Generator, Optional, Sequence, overload
from itertools import islice

class Sieve:
    data: list[int]
    def __init__(self) -> None: ...
    @overload
    def __getitem__(self, item: int) -> int: ...
    @overload
    def __getitem__(self, item: slice) -> list[int]: ...
    def __getitem__(self, item: int | slice) -> int | list[int]: ...
    def __contains__(self, item: int) -> bool: ...
    def __len__(self) -> int: ...
    def __iter__(self) -> Generator: ...
    def search(self, *args: int) -> tuple[int, ...]: ...
    def extend(self, n: int) -> None: ...
    def range(self, a: int, b: Optional[int] = ...) -> Optional[islice]: ...
    @property
    def list(self) -> list[int]: ...
    @property
    def tail(self) -> int: ...

primesieve: Sieve

def miller_rabin(n: int, k: int = ...) -> bool: ...
def miller_rabin_bases(bases: list[int], n: int) -> bool: ...
def baillie_psw(n: int, mr: bool = ...) -> bool: ...
def known_prime(n: int) -> bool: ...
def isprime(n: int) -> bool: ...
def randprime(a: int, b: int = ...) -> int: ...
def confirm_prime(n: int) -> bool: ...
def next_prime(n: int) -> int: ...
def prev_prime(n: int) -> int: ...
def prime_range(a: int, b: Optional[int] = ...) -> list[int]: ...
def sqrt_mod(a: int, p: int) -> Optional[int]: ...
def lift_sqrt(root: int, n: int, modulus: int, q: int = ...) -> int: ...
def quadratic_residue(a: int, p: int) -> bool: ...
def quadratic_non_residue(a: int, p: int) -> bool: ...
def chinese_remainder(values: Sequence[int], moduli: Sequence[int]) -> int: ...
