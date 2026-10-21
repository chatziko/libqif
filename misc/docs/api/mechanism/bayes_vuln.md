# Mechanism construction for Bayes vulnerability

```{py:module} qif.mechanism.bayes_vuln
```

`qif.mechanism.bayes_vuln`

## Functions

```{eval-rst}
.. autofunction:: min_loss_given_max_vuln

   The mechanism with the smallest expected loss among those with posterior Bayes vulnerability
   :math:`V_1[\pi\triangleright C]` at most ``max_vuln``.
```

```{eval-rst}
.. autofunction:: min_vuln_given_max_loss

   The mechanism with the smallest posterior Bayes vulnerability :math:`V_1[\pi\triangleright C]` among those with
   expected loss at most ``max_loss``.
```

```{eval-rst}
.. autofunction:: min_vuln_for_row

   Optimal choice of a single row [ACP25]: ``C`` gives the (fixed) rows of all secrets but one, and ``pi`` their prior
   probabilities (summing to :math:`1-p`). Returns the row :math:`q` for the remaining secret, of prior probability
   ``p``, that minimizes the posterior Bayes vulnerability of the channel obtained by adding :math:`q` to ``C``.
```

## References

[ACP25]
: A. Athanasiou, K. Chatzikokolakis, C. Palamidessi. *Self-Defense: Optimal QIF Solutions and Application to
  Website Fingerprinting*. CSF 2025.
