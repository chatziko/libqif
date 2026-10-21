# Basics

```{py:module} qif
```

`qif`

## Types

Channels and distributions are numpy arrays: a channel $C$ is a matrix with one row per secret and one column
per observation, a prior $\pi$ is a vector. Most functions work both with floating point (`double`) and with
exact rational (`rat`) arithmetic, depending on the type of their arguments.

`double`
: Alias of `numpy.float64`.

`rat`
: Alias of `fractions.Fraction`.

`uint`
: Alias of `numpy.uint32`.

## Functions

```{eval-rst}
.. autofunction:: set_default_type
```

## Classes

```{eval-rst}
.. autoclass:: point
   :members:
```
