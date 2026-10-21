# Shannon entropy

```{py:module} qif.measure.shannon
```

`qif.measure.shannon`

Shannon entropy $H(\pi) = -\sum_x \pi_x \log_2 \pi_x$, an uncertainty measure (in bits). Only `double` is supported.

## Functions

```{eval-rst}
.. autofunction:: prior

   The Shannon entropy :math:`H(\pi)`.
```

```{eval-rst}
.. autofunction:: posterior

   The posterior Shannon entropy :math:`H[\pi\triangleright C]`, that is the conditional entropy :math:`H(X|Y)`.
```

```{eval-rst}
.. autofunction:: add_leakage

   The additive leakage :math:`H(\pi) - H[\pi\triangleright C]`, that is the mutual information :math:`I(X;Y)`.
```

```{eval-rst}
.. autofunction:: mult_leakage

   The multiplicative leakage :math:`H(\pi) / H[\pi\triangleright C]`.
```

```{eval-rst}
.. autofunction:: add_capacity

   The Shannon capacity :math:`\max_\pi I(X;Y)`, computed with the Blahut–Arimoto algorithm. Iterates until the lower
   and upper bounds of the algorithm are equal up to ``md``, ``mrd``. Returns ``(capacity, pi)``, where ``pi`` is a
   prior achieving it.
```
