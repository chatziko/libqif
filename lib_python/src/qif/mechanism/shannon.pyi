"""
Mechanism construction for Shannon enropy.
"""
from .. import typing as t

def max_entropy_given_same_loss(pi: t.ndarray, out: t.ndarray, loss: t.Metric[int,float], md: float = ..., mrd: float = ...) -> t.ndarray: ...
