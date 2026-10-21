# d-privacy

```{py:module} qif.measure.d_privacy
```

`qif.measure.d_privacy`

For a metric $d$ on secrets, a channel $C$ satisfies $d$-privacy if
$C_{x,y} \le e^{d(x,x')} C_{x',y}$ for all secrets $x,x'$ and observations $y$. When the secrets are databases
and $d = \varepsilon d_H$, with $d_H$ the Hamming distance, this is standard $\varepsilon$-differential privacy.
Only `double` is supported.

## Functions

```{eval-rst}
.. autofunction:: prior

   A prior vulnerability, :math:`\max_{x,x'} |\ln \pi_x - \ln \pi_{x'}| / d(x,x')`. Its additive max-case capacity
   (the maximum over all priors of :math:`\max_y V(p_{X|y}) - V(\pi)`) equals :py:func:`smallest_epsilon`.
```

```{eval-rst}
.. autofunction:: is_private

   Whether ``C`` satisfies ``d``-privacy.
```

```{eval-rst}
.. autofunction:: smallest_epsilon

   The smallest :math:`\varepsilon` such that ``C`` satisfies :math:`\varepsilon d`-privacy.
```
