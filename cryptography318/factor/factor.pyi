from collections.abc import Callable
from typing import Optional, Any



def factor(
        n: int, rho: bool = ..., ecm: bool = ..., p1: bool = ..., qs: bool = ..., limit: Optional[int] = ..., *,
        details: dict | None = ...
) -> Optional[dict[int, int]]: ...


def pollard_p1(n: int, B: Optional[int] = ..., _retry: int = ...) -> Optional[int]: ...


def pollard_rho_factor(n: int, mix: Optional[Callable[[int], int]] = ..., _retry: int = ...) -> Optional[int]: ...


def factor_small(factors: dict[int, int], n: int, limit: int) -> tuple[int, int]: ...


def _factor_further(n: int, f: Optional[int], factors: dict[int, int], kwargs) -> bool: ...
