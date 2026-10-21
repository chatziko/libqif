# Mechanism construction for Bayes risk

```{py:module} qif.mechanism.bayes_risk
```

`qif.mechanism.bayes_risk`

## Functions

```{eval-rst}
.. autofunction:: min_loss_given_min_risk

   The mechanism with the smallest expected loss among those with posterior Bayes risk
   :math:`1 - V_1[\pi\triangleright C]` at least ``min_risk``.
```

```{eval-rst}
.. autofunction:: max_risk_given_max_loss

   The mechanism with the largest posterior Bayes risk among those with expected loss at most ``max_loss``.
```

```{eval-rst}
.. autofunction:: max_risk_for_row

   Same as :py:func:`qif.mechanism.bayes_vuln.min_vuln_for_row`: the row that maximizes the posterior Bayes risk is
   the one that minimizes the posterior Bayes vulnerability.
```
