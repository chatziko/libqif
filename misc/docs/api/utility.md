# Utility measures

```{py:module} qif.utility
```

`qif.utility`

The utility of a mechanism $C$ is measured by a distance (or, more generally, a loss function) $d(x,y)$ between the
secret $x$ and the reported observation $y$ (see {doc}`mechanism/index`).

## Functions

```{eval-rst}
.. autofunction:: expected_distance

   The expected distance :math:`\sum_{x,y} \pi_x C_{x,y}\, d(x,y)` between the secret and the reported
   observation. The distance is given either as a matrix ``D`` with :math:`D_{x,y} = d(x,y)`, or as a function
   ``d(x, y)``.
```
