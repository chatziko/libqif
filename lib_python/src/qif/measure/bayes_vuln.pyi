"""
Bayes vulnerability.
"""
from .. import typing as t

def add_leakage(pi: t.ndarray, C: t.ndarray) -> t.FloatOrRat: ...

def min_entropy_leakage(pi: t.ndarray, C: t.ndarray) -> float: ...

def mult_capacity(C: t.ndarray) -> t.FloatOrRat: ...

def mult_leakage(pi: t.ndarray, C: t.ndarray) -> t.FloatOrRat: ...

def posterior(pi: t.ndarray, C: t.ndarray) -> t.FloatOrRat: ...

def prior(pi: t.ndarray) -> t.FloatOrRat: ...

def strategy(pi: t.ndarray, C: t.ndarray) -> t.ndarray: ...

