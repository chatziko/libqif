"""
Mechanism construction for Bayes vulnerability.
"""
from .. import typing as t

def min_loss_given_max_vuln(pi: t.ndarray, n_cols: int, max_vuln: t.FloatOrRat, loss: t.Metric[int,t.FloatOrRat], hard_max_loss: t.FloatOrRat = t.inf) -> t.ndarray: ...

def min_vuln_for_row(pi: t.ndarray, p: t.FloatOrRat, C: t.ndarray) -> t.ndarray: ...

def min_vuln_given_max_loss(pi: t.ndarray, n_cols: int, max_loss: t.FloatOrRat, loss: t.Metric[int,t.FloatOrRat], hard_max_loss: t.FloatOrRat = t.inf) -> t.ndarray: ...

