"""
Utility measures.
"""
from . import typing as t

@t.overload
def expected_distance(D: t.ndarray, pi: t.ndarray, C: t.ndarray) -> t.FloatOrRat: ...
@t.overload
def expected_distance(d: t.Metric[int, t.FloatOrRat], pi: t.ndarray, C: t.ndarray) -> t.FloatOrRat: ...
