"""
Mechanism construction for :math:`g`-vulnerabiliy.
"""
from .. import typing as t

def min_loss_given_max_vuln(pi: t.ndarray, n_cols: int, n_guesses: int, max_vuln: t.FloatOrRat, gain: t.Metric[int,t.FloatOrRat], loss: t.Metric[int,t.FloatOrRat], hard_max_loss: t.FloatOrRat = t.inf) -> t.ndarray: ...

def min_vuln_given_max_loss(pi: t.ndarray, n_cols: int, n_guesses: int, max_loss: t.FloatOrRat, gain: t.Metric[int,t.FloatOrRat], loss: t.Metric[int,t.FloatOrRat], hard_max_loss: t.FloatOrRat = t.inf) -> t.ndarray: ...
