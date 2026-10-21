# Mechanism construction for g-vulnerability

```{py:module} qif.mechanism.g_vuln
```

`qif.mechanism.g_vuln`

The adversary is given by a gain function `gain(w, x)` on `n_guesses` actions $w$ (see {doc}`../measure/g_vuln`).

## Functions

```{eval-rst}
.. autofunction:: min_loss_given_max_vuln

   The mechanism with the smallest expected loss among those with posterior :math:`g`-vulnerability
   :math:`V_g[\pi\triangleright C]` at most ``max_vuln``.
```

```{eval-rst}
.. autofunction:: min_vuln_given_max_loss

   The mechanism with the smallest posterior :math:`g`-vulnerability :math:`V_g[\pi\triangleright C]` among those
   with expected loss at most ``max_loss``.
```
