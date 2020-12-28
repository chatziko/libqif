"""
Refinement relations.
"""
from . import typing as t

# __all__ = [
#     "add_metric",
#     "add_metric_bound",
#     "max_refined_by",
#     "priv_refined_by",
#     "refined_by"
# ]


# @t.overload
def add_metric(pi: t.ndarray, A: t.ndarray, B: t.ndarray) -> t.Tuple[t.rat, t.ndarray]: ...
# @t.overload
# def add_metric(pi: t.ndarray, A: t.ndarray, B: t.ndarray) -> t.Tuple[float, t.ndarray]: ...

# @t.overload
def add_metric_bound(pi: t.ndarray, A: t.ndarray, B: t.ndarray) -> float: ...
# @t.overload
# def add_metric_bound(pi: t.ndarray, A: t.ndarray, B: t.ndarray) -> t.rat: ...

def max_refined_by(A: t.ndarray, B: t.ndarray) -> bool: ...

def priv_refined_by(A: t.ndarray, B: t.ndarray) -> bool: ...

def refined_by(A: t.ndarray, B: t.ndarray, method: str = 'factorize') -> object: ...

