# Refinement relations

```{py:module} qif.refinement
```

`qif.refinement`

A channel $B$ refines $A$, written $A \sqsubseteq B$, if $B = AR$ for some channel $R$, that is if $B$ can be
obtained by post-processing the output of $A$. Equivalently, $B$ never leaks more than $A$:
$V_g[\pi\triangleright A] \ge V_g[\pi\triangleright B]$ for all priors $\pi$ and gain functions $g$. This module
checks refinement and some variants of it, and measures how much $B$ can leak more than $A$ when refinement fails.

## Functions

```{eval-rst}
.. autofunction:: refined_by

   Whether :math:`A \sqsubseteq B`, by checking whether :math:`B = AR` for some channel :math:`R` using
   :py:func:`qif.channel.factorize` with the given ``method``.

   ``method="factorize"`` (same as ``"auto"``) and ``method="project"`` (same as :py:func:`refined_by_witness`
   with ``"project_qp"``) are deprecated.
```

```{eval-rst}
.. autofunction:: refined_by_witness

   Whether :math:`A \sqsubseteq B`, also giving a witness for the answer. Returns ``(refined, G, R)``: if ``refined``
   is true, then :math:`AR = B` (up to the precision of the method) and ``G`` is empty. Otherwise ``R`` is empty and
   ``G`` is a counter-example: a gain matrix (with values in :math:`[0,1]`) such that
   :math:`V_G[\vartheta\triangleright A] < V_G[\vartheta\triangleright B]` for the uniform prior :math:`\vartheta`.
   The counter-example is checked directly, so a false answer is always correct (up to rounding errors).

   Only ``method="project_qp"`` (the default ``"auto"``) gives a witness, so this is ``double`` only.
```

```{eval-rst}
.. autofunction:: max_refined_by

   Whether :math:`A` is refined by :math:`B` in the max-case [CFP19]: :math:`R A^* = B^*` for some channel :math:`R`,
   where :math:`C^*` is obtained from :math:`C` by normalizing its columns and transposing (its rows are the
   posteriors of :math:`C` under the uniform prior). Equivalently, :math:`\max_y V(p_{X|y})` is never larger for
   :math:`B` than for :math:`A`, for all priors and all continuous quasi-convex vulnerabilities :math:`V`.
   ``method``: see :py:func:`qif.channel.factorize`.
```

```{eval-rst}
.. autofunction:: priv_refined_by

   Whether :math:`A` is refined by :math:`B` with respect to privacy [CFP19]: for every metric :math:`d`, if :math:`A`
   satisfies :math:`d`-privacy then so does :math:`B`. Checked by testing whether :math:`B` satisfies
   :math:`d_A`-privacy, where :math:`d_A(x,x')` is the multiplicative distance between the rows :math:`A_{x,-}` and
   :math:`A_{x',-}` (see :doc:`measure/d_privacy`). Only ``double`` is supported.
```

```{eval-rst}
.. autofunction:: add_metric

   The additive refinement metric [CS19]: the maximum of :math:`V_g[\pi\triangleright B] - V_g[\pi\triangleright A]`
   over all one-bounded gain functions :math:`g` (:math:`g \le 1`), that is how much the additive leakage of :math:`B`
   can exceed that of :math:`A` on the prior ``pi``. It is :math:`0` when :math:`A \sqsubseteq B`. Computed by linear
   programming; returns ``(value, G)``, where ``G`` is a gain matrix achieving the maximum.
```

```{eval-rst}
.. autofunction:: add_metric_bound

   An upper bound of :py:func:`add_metric`, cheaper to compute: the Kantorovich distance between the hypers
   :math:`[\pi\triangleright A]` and :math:`[\pi\triangleright B]`, with respect to the quasiconvex separation
   :math:`q` between their inners (see :py:func:`qif.metric.convex_separation_quasi`) [CS19].
```

## References

[CFP19]
: K. Chatzikokolakis, N. Fernandes, C. Palamidessi. *Comparing systems: max-case refinement orders and
  application to differential privacy*. CSF 2019.

[CS19]
: K. Chatzikokolakis, G. Smith. *Refinement metrics for quantitative information flow*. In The Art of Modelling
  Computational Systems, LNCS 11760, 2019.
