# Probability distributions

```{py:module} qif.probab
```

`qif.probab`

A distribution $\pi$ on a finite set $\mathcal X = \{0,\ldots,n{-}1\}$ is a 1-dimensional numpy array of
length $n$, with `pi[x]` $= \pi_x$. Elements are `double` or `rat`; the functions that construct
distributions take a `type` argument, defaulting to the one given to {py:func}`qif.set_default_type`.

For `double` elements, comparisons are made up to a maximum relative difference `mrd` (default `def_mrd`
$= 100\epsilon$, with $\epsilon$ the machine epsilon); `mrd = 0` means exact comparison. For `rat` elements the
default is `0`.

## Functions

```{eval-rst}
.. autofunction:: uniform

   The uniform distribution :math:`\vartheta` on ``n_elem`` elements, that is :math:`\vartheta_x = 1/n`.
```

```{eval-rst}
.. autofunction:: point

   The point distribution :math:`[x]`, giving probability :math:`1` to ``x``.
```

```{eval-rst}
.. autofunction:: randu

   A distribution chosen uniformly at random from :math:`\mathbb D\mathcal X` (the probability simplex).
```

```{eval-rst}
.. autofunction:: normalize

   Returns ``pi`` scaled to sum to :math:`1`, that is :math:`\pi_x / \sum_{x'} \pi_{x'}`.
```

```{eval-rst}
.. autofunction:: sample

   Draws an element :math:`x` from ``pi``. Returns a single element,
   or an array of ``n_samples`` independent samples.
```

```{eval-rst}
.. autofunction:: is_uniform

   Whether ``pi`` is the uniform distribution :math:`\vartheta`, up to ``mrd``.
```

```{eval-rst}
.. autofunction:: is_proper

   Whether ``pi`` is a proper distribution: its elements are non-negative and sum to :math:`1` (up to ``mrd``).
```

```{eval-rst}
.. autofunction:: assert_proper

   Raises ``RuntimeError`` if ``pi`` is not a proper distribution (see :py:func:`is_proper`).
```

```{eval-rst}
.. autofunction:: to_grid

   For a distribution ``pi`` on the cells of a grid with ``width`` columns, where cell :math:`(x,y)` has index
   ``y*width + x`` (as in :py:meth:`qif.point.from_cell`), returns the grid as a matrix ``M`` with ``M[y, x]`` being
   the probability of cell :math:`(x,y)`.
```

```{eval-rst}
.. autofunction:: from_grid

   The inverse of :py:func:`to_grid`: returns the distribution ``pi`` with ``pi[y*width + x] = grid[y, x]``, where
   ``width`` is the number of columns of ``grid``.
```
