"""
Mechanism construction for :math:`\ell`-risk.
"""
from .. import typing as t

def max_risk_given_max_loss(pi: t.ndarray, n_cols: int, n_guesses: int, max_loss: t.FloatOrRat, adv_loss: t.Metric[int,t.FloatOrRat], loss: t.Metric[int,t.FloatOrRat], hard_max_loss: t.FloatOrRat = t.inf) -> t.ndarray: ...

def min_loss_given_min_risk(pi: t.ndarray, n_cols: int, n_guesses: int, min_risk: t.FloatOrRat, adv_loss: t.Metric[int,t.FloatOrRat], loss: t.Metric[int,t.FloatOrRat], hard_max_loss: t.FloatOrRat = t.inf) -> t.ndarray: ...
