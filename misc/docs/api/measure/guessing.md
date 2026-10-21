# Guessing entropy

```{py:module} qif.measure.guessing
```

`qif.measure.guessing`

Guessing entropy $G(\pi) = \sum_{i=1}^n i\, \pi_{x_i}$, where $x_1, \ldots, x_n$ are the secrets in non-increasing
order of probability: the expected number of guesses the adversary needs to find the secret, guessing one value
at a time in the best order. It is an uncertainty measure.

## Functions

```{eval-rst}
.. autofunction:: prior

   The guessing entropy :math:`G(\pi)`.
```

```{eval-rst}
.. autofunction:: posterior

   The posterior guessing entropy :math:`G[\pi\triangleright C] = \sum_y p(y)\, G(p_{X|y})`.
```

```{eval-rst}
.. autofunction:: add_leakage

   The additive leakage :math:`G(\pi) - G[\pi\triangleright C]`.
```

```{eval-rst}
.. autofunction:: mult_leakage

   The multiplicative leakage :math:`G(\pi) / G[\pi\triangleright C]`.
```
