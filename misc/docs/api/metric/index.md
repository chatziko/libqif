# Metrics and other distance functions

```{py:module} qif.metric
```

`qif.metric`

A metric is represented as a function `d(a, b)` returning a number; any Python function can be used. The functions
of this module construct common metrics, and combine existing ones into new ones. Metrics can be defined on
secrets (indices $0,\ldots,n{-}1$), on numbers, on {py:class}`qif.point`, or on distributions. Some of the
constructed functions are not strictly metrics (for instance they are not symmetric), but are still useful as
distances.

```{toctree}
:maxdepth: 1

optimize
```

## Functions

```{eval-rst}
.. autofunction:: euclidean

   The Euclidean metric :math:`|a - b|` on numbers (``double``, ``rat``) or secrets (``uint``), or the Euclidean
   distance between two points (``point``), depending on ``type``.
```

```{eval-rst}
.. autofunction:: euclidean_chain

   The "chainable" predicate of :py:func:`euclidean` for the given ``type``: ``chain(a, b)`` is true if there is a
   chain :math:`a = e_1, \ldots, e_k = b` of intermediate elements with :math:`d(a,b) = \sum_i d(e_i, e_{i+1})`. For
   secrets (``uint``) this holds for all :math:`a, b` with :math:`|a - b| \neq 1`; for other types it is always false.
   Used by :py:func:`qif.mechanism.d_privacy.min_loss_given_d` to skip redundant constraints.
```

```{eval-rst}
.. autofunction:: discrete

   The discrete metric: :math:`0` if :math:`a = b`, :math:`1` otherwise.
```

```{eval-rst}
.. autofunction:: mult_reals

   The metric :math:`|\ln a - \ln b|` on positive numbers.
```

```{eval-rst}
.. autofunction:: scale

   The metric :math:`\mathrm{coeff} \cdot d`. Distances equal to :math:`0` stay :math:`0` even when ``coeff`` is infinite.
```

```{eval-rst}
.. autofunction:: min

   The pointwise minimum :math:`\min(d_1(a,b), d_2(a,b))` (not necessarily a metric).
```

```{eval-rst}
.. autofunction:: max

   The pointwise maximum :math:`\max(d_1(a,b), d_2(a,b))`.
```

```{eval-rst}
.. autofunction:: mirror

   The distance :math:`d(b,a)`, with the arguments swapped (useful for non-symmetric distances).
```

```{eval-rst}
.. autofunction:: threshold

   The distance :math:`\min(d(a,b), \mathrm{thres})`.
```

```{eval-rst}
.. autofunction:: threshold_inf

   The distance :math:`d(a,b)` if it is at most ``thres``, infinity otherwise.
```

```{eval-rst}
.. autofunction:: from_distance_matrix

   The distance on secrets given by a matrix: :math:`d(i,j) = M_{i,j}`.
```

```{eval-rst}
.. autofunction:: to_distance_matrix

   The ``n_rows`` :math:`\times` ``n_cols`` matrix :math:`M_{i,j} = d(i,j)` of a distance on secrets (``n_cols = 0``
   means ``n_cols = n_rows``). Applied to a gain function :math:`g(w,x)`, it gives the gain matrix.
```

```{eval-rst}
.. autofunction:: l1

   The :math:`\ell_1` metric :math:`\sum_i |a_i - b_i|` on vectors.
```

```{eval-rst}
.. autofunction:: l2

   The Euclidean (:math:`\ell_2`) metric :math:`\sqrt{\sum_i (a_i - b_i)^2}` on vectors.
```

```{eval-rst}
.. autofunction:: total_variation

   The total variation distance :math:`\frac12 \sum_x |a_x - b_x|` between distributions, equal to the largest
   difference :math:`|a(Z) - b(Z)|` in the probability of an event :math:`Z`.
```

```{eval-rst}
.. autofunction:: mult_total_variation

   The multiplicative distance :math:`\max_x |\ln a_x - \ln b_x|` between distributions (infinity if their supports
   differ), equal to the largest :math:`|\ln (a(Z)/b(Z))|` for an event :math:`Z`. A channel :math:`C` satisfies
   :math:`d`-privacy iff the multiplicative distance between rows :math:`C_{x,-}, C_{x',-}` is at most
   :math:`d(x,x')` (see :doc:`../measure/d_privacy`).
```

```{eval-rst}
.. autofunction:: convex_separation_quasi

   The quasiconvex separation :math:`q(a,b) = \max_{x \in \lceil a\rceil} (1 - b_x / a_x)`, where :math:`\lceil a\rceil` is
   the support of :math:`a` (a non-symmetric distance). The additive capacity for a fixed prior :math:`\pi` over
   one-bounded gain functions is the expected value of :math:`q(\sigma, \pi)` over the posteriors :math:`\sigma` of
   :math:`[\pi\triangleright C]` (see :py:func:`qif.measure.g_vuln.add_capacity`).
```

```{eval-rst}
.. autofunction:: convex_separation

   The symmetric version :math:`\max(q(a,b), q(b,a))` of :py:func:`convex_separation_quasi`.
```

```{eval-rst}
.. autofunction:: kantorovich

   The Kantorovich lifting of a metric ``d`` on secrets to a metric on distributions: the minimum cost of transforming
   :math:`a` into :math:`b` by moving probability mass, where moving a unit of mass from :math:`i` to :math:`j` costs
   :math:`d(i,j)`. Equivalently, the maximum of :math:`|\sum_x a_x f(x) - \sum_x b_x f(x)|` over all functions
   :math:`f` with :math:`|f(x) - f(x')| \le d(x,x')`. With ``method="fastemd"`` (the default ``"auto"``) it is computed
   with the FastEMD algorithm, which requires ``d`` to be a metric and has limited precision (about ``1e-5``
   relative). ``method="lp"`` uses linear programming instead: slower, but more precise, and ``d`` can be any
   distance function.
```

```{eval-rst}
.. autofunction:: mult_kantorovich

   The multiplicative variant of the Kantorovich lifting [CGPX14]: the maximum of
   :math:`|\ln \sum_x a_x f(x) - \ln \sum_x b_x f(x)|` over all functions :math:`f` with values in :math:`[0,1]` and
   :math:`|\ln f(x) - \ln f(x')| \le d(x,x')`. Computed by linear programming (via the dual form).
```

## References

[CGPX14]
: K. Chatzikokolakis, D. Gebler, C. Palamidessi, L. Xu. *Generalized bisimulation metrics*. CONCUR 2014.
