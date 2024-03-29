from typing import Optional


def baby_step_giant_step(g: int, h: int, p: int, order: Optional[int] = ...) -> Optional[int]: ...


def pollard_rho_dlp(g: int, h: int, p: int, order: Optional[int] = ...) -> int: ...


def calculate_state(state: tuple[int, int, int], g: int, h: int, p: int) -> tuple[int, int, int]: ...


def pohlig_hellman(g: int, h: int, p: int, q: int, exp, prog: bool = ...) -> int: ...
