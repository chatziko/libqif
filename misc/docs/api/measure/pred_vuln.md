# Vulnerability of guessing a predicate of the secret

```{py:module} qif.measure.pred_vuln
```

`qif.measure.pred_vuln`

The vulnerability of the adversary guessing whether the secret satisfies a predicate $P$, given as a $0/1$ array `P`
over the secrets. This is $g$-vulnerability for the two-block gain function, whose two actions are "$x \in P$" and
"$x \notin P$" (see {doc}`g_vuln`); for instance $V(\pi) = \max(\pi(P), \pi(\overline P))$.

## Functions

```{eval-rst}
.. autofunction:: G_pred

   The gain matrix of the two-block gain function: row :math:`0` is the action ":math:`x \in P`" (gain :math:`1` if the
   secret satisfies :math:`P`), row :math:`1` is ":math:`x \notin P`".
```

```{eval-rst}
.. autofunction:: prior

   The prior vulnerability :math:`\max(\pi(P), \pi(\overline P))`.
```

```{eval-rst}
.. autofunction:: posterior

   The posterior vulnerability.
```

```{eval-rst}
.. autofunction:: add_leakage

   The additive leakage.
```

```{eval-rst}
.. autofunction:: mult_leakage

   The multiplicative leakage.
```

```{eval-rst}
.. autofunction:: mult_capacity

   The multiplicative capacity: the maximum of :py:func:`mult_leakage` over all priors. Returns ``(capacity, pi)``,
   where ``pi`` is a prior achieving it.

   The capacity is :math:`1 + \frac12 \max_{x \in P, x' \notin P} \lVert C_{x,-} - C_{x',-}\rVert_1`, achieved by the
   uniform prior on :math:`x, x'`. Only ``method="direct"`` is available.
```

```{eval-rst}
.. autofunction:: binary_channel

   Reduces the problem to a binary secret ":math:`P`" / ":math:`\overline P`". Returns ``(pi_bin, C_bin)``, where
   ``pi_bin`` :math:`= (\pi(P), \pi(\overline P))` and the two rows of ``C_bin`` are the average of the rows of
   ``C`` on :math:`P` and on :math:`\overline P` (weighted by ``pi``). The vulnerability of :math:`P` under ``(pi, C)``
   equals the Bayes vulnerability under ``(pi_bin, C_bin)``.
```
