"""
Channel composition.
"""
from .. import typing as t

def parallel(A: t.ndarray, B: t.ndarray) -> t.ndarray: ...

def repeated_independent(C: t.ndarray, n: int) -> t.ndarray: ...
