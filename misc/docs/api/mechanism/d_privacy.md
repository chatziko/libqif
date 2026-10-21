# Mechanism construction for d-privacy

```{py:module} qif.mechanism.d_privacy
```

`qif.mechanism.d_privacy`

Mechanisms satisfying $d$-privacy (see {doc}`../measure/d_privacy`). The privacy parameter $\varepsilon$ is
either given explicitly, or included in the metric (for instance by using {py:func}`qif.metric.scale` to get
$\varepsilon d$). Only `double` is supported.

## Functions

```{eval-rst}
.. autofunction:: distance_matrix

   The ``n_rows`` :math:`\times` ``n_cols`` matrix :math:`e^{-d(i,j)}`, the building block of most mechanisms in this
   module.
```

```{eval-rst}
.. autofunction:: geometric

   The truncated geometric mechanism [GRS09], satisfying :math:`\varepsilon d`-privacy for :math:`d(x,x') = |x-x'|`.
   Row :math:`i` corresponds to the value ``first_x + i`` and column :math:`j` to the value ``first_y + j``
   (``n_cols = 0`` means ``n_cols = n_rows``). The probability of reporting :math:`y` for :math:`x` is proportional to
   :math:`e^{-\varepsilon|x-y|}`, and the probability of values outside the range of columns is moved to the first or
   last column.
```

```{eval-rst}
.. autofunction:: exponential

   The exponential mechanism: :math:`C_{x,y}` is proportional to :math:`e^{-d(x,y)/2}`, where :math:`d` is a metric
   between secrets and observations (``n_cols = 0`` means ``n_cols = n_rows``). It satisfies :math:`d`-privacy.
```

```{eval-rst}
.. autofunction:: randomized_response

   The randomized response mechanism: :math:`C_{x,y}` is proportional to :math:`1` if :math:`y = x` and to
   :math:`e^{-\varepsilon}` otherwise (``n_cols = 0`` means ``n_cols = n_rows``). It satisfies
   :math:`\varepsilon`-differential privacy between any two secrets.
```

```{eval-rst}
.. autofunction:: tight_constraints

   The tight-constraints mechanism for the metric ``d`` [ECP13]: the mechanism with :math:`C_{x,y} = c_y e^{-d(x,y)}`
   for some non-negative coefficients :math:`c_y`, which is optimal (wrt utility) for a whole region of priors.
   Returns an empty array if no such mechanism exists.
```

```{eval-rst}
.. autofunction:: exact_distance

   A mechanism (with ``n_rows`` :math:`+1` columns) such that the multiplicative distance between any two rows
   :math:`C_{x,-}, C_{x',-}` is exactly :math:`d(x,x')` [CFP19]: it satisfies :math:`d`-privacy, but no
   :math:`d'`-privacy for a smaller metric :math:`d'`.
```

```{eval-rst}
.. autofunction:: min_loss_given_d

   The mechanism with the smallest expected ``loss`` among those satisfying ``d_priv``-privacy, computed by linear
   programming [BCP14]. Returns an empty array if the linear program fails.

   ``vars`` controls which probabilities are variables of the program. With ``"all"`` (default) every
   :math:`C_{x,y}` is a variable. With ``"dist"`` (resp. ``"dist_strict"``), :math:`C_{x,y}` is restricted to depend
   only on :math:`d_\mathrm{priv}(x,y)` and :math:`y` (resp. only on :math:`d_\mathrm{priv}(x,y)`), with secrets
   and observations in the same space; this gives a much smaller program, but the result is not guaranteed to be
   optimal.

   With ``"all"``, the remaining parameters reduce the size of the program: ``d_priv_ch`` is a "chainable" predicate
   for ``d_priv`` (see :py:func:`qif.metric.euclidean_chain`), the constraints of chainable pairs being redundant,
   and distances larger than ``inf`` impose no constraint.
```

## References

[BCP14]
: N. E. Bordenabe, K. Chatzikokolakis, C. Palamidessi. *Optimal geo-indistinguishable mechanisms for location
  privacy*. CCS 2014.

[CFP19]
: K. Chatzikokolakis, N. Fernandes, C. Palamidessi. *Comparing systems: max-case refinement orders and
  application to differential privacy*. CSF 2019.

[ECP13]
: E. ElSalamouny, K. Chatzikokolakis, C. Palamidessi. *A differentially private mechanism of optimal utility for a
  region of priors*. POST 2013.

[GRS09]
: A. Ghosh, T. Roughgarden, M. Sundararajan. *Universally utility-maximizing privacy mechanisms*. STOC 2009.
