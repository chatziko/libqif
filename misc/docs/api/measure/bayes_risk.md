# Bayes risk

```{py:module} qif.measure.bayes_risk
```

`qif.measure.bayes_risk`

Bayes risk $1 - V_1(\pi)$ is the probability that the adversary guesses the secret wrongly in one try (see
{doc}`bayes_vuln`). It is an uncertainty measure: lower values are worse for the user.

## Functions

```{eval-rst}
.. autofunction:: prior

   The Bayes risk :math:`1 - V_1(\pi)`.
```

```{eval-rst}
.. autofunction:: posterior

   The posterior Bayes risk :math:`1 - V_1[\pi\triangleright C]`.
```

```{eval-rst}
.. autofunction:: add_leakage

   The additive leakage: prior minus posterior risk (equal to the additive Bayes leakage).
```

```{eval-rst}
.. autofunction:: mult_leakage

   The multiplicative leakage: prior risk divided by posterior risk. It is :math:`1` if the prior risk is
   :math:`0`, and infinity if only the posterior risk is :math:`0`.
```

```{eval-rst}
.. autofunction:: mult_capacity

   The multiplicative capacity: the maximum of :py:func:`mult_leakage` over all priors. Returns ``(capacity, pi)``,
   where ``pi`` is a prior achieving it.

   The capacity is :math:`1/(1 - \frac12 \max_{x,x'} \lVert C_{x,-} - C_{x',-}\rVert_1)`, achieved by the uniform
   prior on the two secrets :math:`x, x'` whose rows are furthest apart. With ``method="direct"`` all pairs of rows are
   compared; ``method="linf"`` uses an embedding in :math:`\ell_\infty`, faster when ``C`` has few columns.
```

```{eval-rst}
.. autofunction:: strategy

   The optimal guess for each observation (same as :py:func:`qif.measure.bayes_vuln.strategy`).
```
