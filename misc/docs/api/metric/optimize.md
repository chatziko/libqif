# Metric optimization problems

```{py:module} qif.metric.optimize
```

`qif.metric.optimize`

## Functions

```{eval-rst}
.. autofunction:: l1_diameter

   The :math:`\ell_1` diameter of the rows of ``C``: the maximum :math:`\lVert C_{x,-} - C_{x',-}\rVert_1`. Returns
   ``(diameter, x, x')`` with the two rows achieving it. With ``method="direct"`` all pairs of rows are compared;
   ``method="linf"`` embeds the rows in :math:`\ell_\infty` of dimension :math:`2^m` (:math:`m` the number of
   columns), faster when ``C`` has few columns.
```

```{eval-rst}
.. autofunction:: l2_min_enclosing_ball

   The smallest Euclidean ball containing the rows of ``C`` [ACP25]. Returns ``(radius, center)``. Only ``double`` is
   supported.
```

```{eval-rst}
.. autofunction:: simplex_l1_min_enclosing_ball

   The distribution :math:`q` minimizing the maximum :math:`\ell_1` distance :math:`\lVert C_{x,-} - q\rVert_1` from the
   rows of ``C`` [ACP25]. Returns ``(radius, q)``. Computed by linear programming, using either ``method="lp"`` or
   ``method="linf"`` (an embedding in :math:`\ell_\infty`, as in :py:func:`l1_diameter`). If ``in_conv_hull`` is true
   (only with ``"lp"``), :math:`q` is restricted to the convex hull of the rows.
```

```{eval-rst}
.. autofunction:: simplex_project

   The Euclidean projection of the vector ``pi`` on the probability simplex, that is the distribution closest to
   ``pi`` in :math:`\ell_2` distance.
```

## References

[ACP25]
: A. Athanasiou, K. Chatzikokolakis, C. Palamidessi. *Self-Defense: Optimal QIF Solutions and Application to
  Website Fingerprinting*. CSF 2025.
