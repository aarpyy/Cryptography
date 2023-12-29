from random import Random
from typing import Self

rndm: Random


TPoint = Identity | WeierstrassPoint


class Identity:
    def __eq__(self, other: TPoint) -> TPoint: ...

    def __neg__(self): ...

    def __add__(self, other: TPoint) -> TPoint: ...

    def __mul__(self, other): ...


class WeierstrassPoint:
    infinity: Identity
    x: int
    y: int
    n: int
    a: int

    def __new__(cls, x, y, a, n): ...

    def __eq__(self, other: TPoint) -> bool: ...

    def __neg__(self) -> TPoint: ...

    def __add__(self, other: TPoint) -> TPoint: ...

    def __mul__(self, other: TPoint) -> TPoint: ...


def gen_weierstrass_point(n: int) -> WeierstrassPoint: ...


class MontgomeryPoint:
    x: int
    z: int
    a_24: int
    n: int

    def __new__(cls, x: int, z: int, a_24: int, n: int) -> MontgomeryPoint: ...

    def add(self, Q: Self, difference: Self) -> Self: ...

    def double(self) -> Self: ...

    def ladder(self, k: int) -> Self: ...


class MontgomeryPointZ1(MontgomeryPoint):
    def __new__(cls, x: int, z: int, a_24: int, n: int) -> MontgomeryPointZ1: ...

    def add(self, Q: Self, difference: Self) -> Self: ...


def make_point(x: int, z: int, a_24: int, n: int) -> MontgomeryPoint: ...


def ecm_weierstrass(N: int, B: int | None = ..., retry: int = ..., verbose: bool = ...) -> int | None: ...


def ecm_mont_phase1(N: int, use_z1: bool = ...) -> int | MontgomeryPoint: ...


def ecm_mont(N, B1: int | None = ..., B2: int | None = ..., retry: int = ..., verbose: bool = ...,
              use_z1: bool = ...) -> int | None: ...


def ecm_mont_basic(N: int, B1: int | None = ..., B2: int | None = ..., retry: int = ...) -> int | None: ...


def ecm(N: int, B1: int | None = ..., B2: int | None = ..., retry: int = ..., verbose: bool = ..., use_z1: bool = ...,
        use_weierstrass: bool = ..., use_basic: bool = ...) -> int | None: ...
