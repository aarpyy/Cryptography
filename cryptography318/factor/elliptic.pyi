from typing import Optional

from .curve import Montgomery, MontgomeryPoint, Weierstrass, WeierstrassPoint


def safe_weierstrass(p: int) -> tuple[Weierstrass, WeierstrassPoint]: ...


def safe_montgomery(p: int) -> tuple[Montgomery, MontgomeryPoint] | int: ...


def z_1_montgomery(p: int) -> tuple[Montgomery, MontgomeryPoint] | int: ...


def ecm_mont_basic(N: int, *, B1: int = ..., B2: int = ..., retry: int = ...) -> Optional[int]: ...


def ecm_mont(N: int, *, B1: int = ..., B2: int = ..., retry: int = ...) -> Optional[int]: ...


def ecm_weierstrass(N: int, *, B: int = ..., retry: int = ...) -> Optional[int]: ...


def lenstra_ecm(N: int, *, B: int = ..., retry: int = ...) -> Optional[int]: ...
