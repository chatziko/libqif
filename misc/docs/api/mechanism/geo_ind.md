# Mechanism construction for geo-indistinguishability

```{py:module} qif.mechanism.geo_ind
```

`qif.mechanism.geo_ind`

Mechanisms for geo-indistinguishability [ABCP13], that is $d$-privacy for $\varepsilon$ times the Euclidean
distance between locations. Locations are {py:class}`qif.point` objects. On a grid of `width` $\times$ `height`
cells with spacing `step`, the cell $(x,y)$ has index `y*width + x` (see {py:meth}`qif.point.from_cell`). Only
`double` is supported.

## Functions

```{eval-rst}
.. autofunction:: planar_laplace_sample

   A random noise vector from the planar Laplace distribution, with density proportional to
   :math:`e^{-\varepsilon \lVert z\rVert}`. Adding it to the real location gives the reported one [ABCP13].
```

```{eval-rst}
.. autofunction:: planar_laplace_grid

   The planar Laplace mechanism on a grid: the probability of reporting each cell, computed by numerical integration.
   Reported locations outside the grid are moved to the closest cell on its boundary.
```

```{eval-rst}
.. autofunction:: planar_geometric_sample

   A random noise vector from the planar geometric distribution, a discrete version of planar Laplace on a grid with
   spacing ``cell_size``: the probability of :math:`z` is proportional to :math:`e^{-\varepsilon \lVert z\rVert}`
   [CEP17]. Returns one vector, or a list of ``n_samples`` independent ones.
```

```{eval-rst}
.. autofunction:: planar_geometric_grid

   The planar geometric mechanism on a grid; reported locations outside the grid are moved to the closest cell on its
   boundary.
```

## References

[ABCP13]
: M. E. Andrés, N. E. Bordenabe, K. Chatzikokolakis, C. Palamidessi. *Geo-indistinguishability: Differential
  privacy for location-based systems*. CCS 2013.

[CEP17]
: K. Chatzikokolakis, E. ElSalamouny, C. Palamidessi. *Efficient utility improvement for location privacy*.
  PoPETs 2017.
