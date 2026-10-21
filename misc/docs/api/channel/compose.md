# Channel composition

```{py:module} qif.channel.compose
```

`qif.channel.compose`

## Functions

```{eval-rst}
.. autofunction:: parallel

   The parallel composition :math:`A \parallel B` of two channels on the same secrets: both are run
   independently on the same secret and both observations are revealed, that is
   :math:`(A\parallel B)_{x,(y_1,y_2)} = A_{x,y_1} B_{x,y_2}`. The observation :math:`(y_1,y_2)` is column
   ``y1 * B.shape[1] + y2``.
```

```{eval-rst}
.. autofunction:: repeated_independent

   The parallel composition of ``n`` copies of ``C``, modelling ``n`` independent runs
   of ``C`` on the same secret (for ``n = 0``, a channel that leaks nothing).
```
