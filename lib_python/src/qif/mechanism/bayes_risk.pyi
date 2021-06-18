"""
Mechanism construction for Bayes risk.
"""
from .. import typing as t

def max_risk_for_row(pi: t.ndarray, p: t.FloatOrRat, C: t.ndarray) -> t.ndarray: ...

def max_risk_given_max_loss(pi: t.ndarray, n_cols: int, max_loss: t.FloatOrRat, loss: t.Metric[int,t.FloatOrRat], hard_max_loss: t.FloatOrRat = t.inf) -> t.ndarray: ...

def min_loss_given_min_risk(pi: t.ndarray, n_cols: int, min_risk: t.FloatOrRat, loss: t.Metric[int,t.FloatOrRat], hard_max_loss: t.FloatOrRat = t.inf) -> t.ndarray: ...

