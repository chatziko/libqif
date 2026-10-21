# Bayes vulnerability

```{py:module} qif.measure.bayes_vuln
```

`qif.measure.bayes_vuln`

Bayes vulnerability $V_1(\pi) = \max_x \pi_x$ is the probability that the adversary guesses the secret
correctly in one try.

## Functions

```{eval-rst}
.. autofunction:: prior

   The Bayes vulnerability :math:`V_1(\pi) = \max_x \pi_x`.
```

```{eval-rst}
.. autofunction:: posterior

   The posterior Bayes vulnerability :math:`V_1[\pi\triangleright C]`, that is the sum of the column maxima of
   the joint matrix :math:`\sum_y \max_x \pi_x C_{x,y}`.
```

```{eval-rst}
.. autofunction:: add_leakage

   The additive Bayes leakage :math:`\mathcal L^+_1(\pi, C) = V_1[\pi\triangleright C] - V_1(\pi)`.
```

```{eval-rst}
.. autofunction:: mult_leakage

   The multiplicative Bayes leakage :math:`\mathcal L^\times_1(\pi, C) = V_1[\pi\triangleright C] / V_1(\pi)`.
```

```{eval-rst}
.. autofunction:: min_entropy_leakage

   The min-entropy leakage :math:`\log_2 \mathcal L^\times_1(\pi, C)` (in bits).
```

```{eval-rst}
.. autofunction:: mult_capacity

   The multiplicative Bayes capacity :math:`\mathcal{ML}^\times_1(\mathbb D, C) = \max_\pi \mathcal L^\times_1(\pi, C)`,
   equal to the sum of the column maxima :math:`\sum_y \max_x C_{x,y}` (always achieved by the uniform prior).
```

```{eval-rst}
.. autofunction:: strategy

   The optimal guess for each observation: an array whose element :math:`y` is the :math:`x` maximizing
   :math:`\pi_x C_{x,y}`.
```
