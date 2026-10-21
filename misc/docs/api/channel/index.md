# Channels

```{py:module} qif.channel
```

`qif.channel`

A channel $C$ from $\mathcal X = \{0,\ldots,n{-}1\}$ to $\mathcal Y = \{0,\ldots,m{-}1\}$ is an $n\times m$ numpy
matrix with $C_{x,y} = p(y|x)$: each row $C_{x,-}$ is the distribution on observations $y$ produced by the secret
$x$. Elements are `double` or `rat`, as for {doc}`distributions <../probab>`.

For a prior $\pi$, the joint matrix $J = \pi\triangleright C$ is given by $J_{x,y} = \pi_x C_{x,y}$; its column
sums give the probability $p(y)$ of each observation, and normalizing a column $y$ with $p(y)>0$ gives the
posterior distribution $p_{X|y}$.

```{toctree}
:maxdepth: 1

compose
```

## Functions

```{eval-rst}
.. autofunction:: normalize

   Returns ``C`` with each row divided by its sum, so that every row becomes a distribution.
```

```{eval-rst}
.. autofunction:: identity

   The identity channel on ``n_rows`` secrets (the identity matrix), which reveals the secret completely.
```

```{eval-rst}
.. autofunction:: no_interference

   A channel that leaks nothing: every row is the point distribution on observation :math:`0`
   (the remaining ``n_cols``\ :math:`{-}1` columns are zero).
```

```{eval-rst}
.. autofunction:: randu

   A channel with each row chosen uniformly at random from the distributions on ``n_cols``
   observations (``n_cols = 0`` means ``n_cols = n_rows``).
```

```{eval-rst}
.. autofunction:: deterministic

   The deterministic channel that maps each secret :math:`x` to the observation ``map(x)``, that is
   :math:`C_{x,\mathrm{map}(x)} = 1`.
```

```{eval-rst}
.. autofunction:: is_proper

   Whether ``C`` is a proper channel: every row is a distribution (non-negative, summing to :math:`1`
   up to ``mrd``).
```

```{eval-rst}
.. autofunction:: assert_proper

   Raises ``RuntimeError`` if ``C`` is not a proper channel (see :py:func:`is_proper`).
```

```{eval-rst}
.. autofunction:: posterior

   The posterior distribution :math:`p_{X|y}` on secrets after observing :math:`y` = ``col``, under prior
   ``pi``. Requires :math:`p(y) > 0`.
```

```{eval-rst}
.. autofunction:: posteriors

   All posterior distributions: a matrix of the same size as ``C`` whose column :math:`y` is
   :math:`p_{X|y}` (undefined for columns with :math:`p(y) = 0`). If ``pi`` is not given the uniform prior is used.
```

```{eval-rst}
.. autofunction:: hyper

   The hyper-distribution :math:`[\pi\triangleright C]`, as a pair ``(outer, inners)``: ``inners`` has one column
   per inner distribution, and ``outer[i]`` is the probability of inner ``i``. Observations producing the same
   posterior are merged into a single inner.
```

```{eval-rst}
.. autofunction:: reduced

   The reduced form of ``C``: zero columns are removed and similar columns (multiples of each other) are
   merged. Two channels produce the same hyper on every prior iff their reduced forms are equal (up to the order
   of columns).
```

```{eval-rst}
.. autofunction:: iterative_bayesian_update

   Estimates the prior from an observed output distribution ``out``, using the iterative
   Bayesian update (an instance of expectation–maximization). Starts from ``start`` (uniform if not given) and stops
   when the :math:`\ell_1` distance between consecutive estimates is at most ``max_diff``, or after ``max_iter``
   iterations (``0`` means no limit). Returns ``(pi, iterations)``.
```

```{eval-rst}
.. autofunction:: factorize

   Returns a channel :math:`X` such that :math:`A = BX`, or an empty matrix if none exists. If
   ``col_stoch`` is true, :math:`X` is required to be column-stochastic instead. Solved by linear programming
   (or a subgradient method for large ``double`` matrices).
```

```{eval-rst}
.. autofunction:: left_factorize

   Returns a channel :math:`X` such that :math:`A = XB`, or an empty matrix if none exists
   (see :py:func:`factorize`).
```

```{eval-rst}
.. autofunction:: sum_column_min

   The sum of the column minima :math:`\sum_y \min_x C_{x,y}`.
```

```{eval-rst}
.. autofunction:: sample

   Draws a secret :math:`x` from ``pi`` and then an observation :math:`y` from :math:`C_{x,-}`. Returns
   the pair ``(x, y)``, or an ``n_samples``\ :math:`\times 2` array of independent samples.
```
