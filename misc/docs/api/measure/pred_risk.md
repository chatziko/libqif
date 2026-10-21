# Risk of guessing a predicate of the secret

```{py:module} qif.measure.pred_risk
```

`qif.measure.pred_risk`

The risk of the adversary guessing wrongly whether the secret satisfies a predicate $P$, given as a $0/1$ array
`P` over the secrets. This is $\ell$-risk for the two-block loss function (see {doc}`l_risk`), equal to
$1$ minus the vulnerability of {doc}`pred_vuln`; for instance $U(\pi) = \min(\pi(P), \pi(\overline P))$.

## Functions

```{eval-rst}
.. autofunction:: L_pred

   The loss matrix of the two-block loss function: row :math:`0` is the action ":math:`x \notin P`" (loss :math:`1` if
   the secret satisfies :math:`P`), row :math:`1` is ":math:`x \in P`".
```

```{eval-rst}
.. autofunction:: prior

   The prior risk :math:`\min(\pi(P), \pi(\overline P))`.
```

```{eval-rst}
.. autofunction:: posterior

   The posterior risk.
```

```{eval-rst}
.. autofunction:: add_leakage

   The additive leakage (prior minus posterior risk).
```

```{eval-rst}
.. autofunction:: mult_leakage

   The multiplicative leakage (prior divided by posterior risk).
```

```{eval-rst}
.. autofunction:: mult_capacity

   The multiplicative capacity: the maximum of :py:func:`mult_leakage` over all priors. Returns ``(capacity, pi)``,
   where ``pi`` is a prior achieving it.

   The capacity is :math:`1/(1 - \frac12 \max_{x \in P, x' \notin P} \lVert C_{x,-} - C_{x',-}\rVert_1)`, achieved by
   the uniform prior on :math:`x, x'`. Only ``method="direct"`` is available.
```

```{eval-rst}
.. autofunction:: binary_channel

   Same as :py:func:`qif.measure.pred_vuln.binary_channel`.
```
