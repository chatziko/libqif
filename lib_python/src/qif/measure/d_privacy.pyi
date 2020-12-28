"""
:math:`d`-privacy.
"""
from .. import typing as t

def is_private(C: t.ndarray, d: t.Metric[int,float]) -> bool: ...

def prior(pi: t.ndarray, d: t.Metric[int,float]) -> float: ...

def smallest_epsilon(C: t.ndarray, d: t.Metric[int,float]) -> float: ...

