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
: Alias of `numpy.uint32`, used for secrets (indices) in metrics, for instance `metric.euclidean(uint)`.

## Functions

```{eval-rst}
.. autofunction:: set_default_type

   Sets the default element type (``double`` or ``rat``) for the functions that construct distributions, channels
   or gain functions, used when their ``type`` argument is not given. Initially ``double``.
```

## Classes

```{eval-rst}
.. py:class:: point(x, y)

   A point in the plane, used for locations (see :doc:`mechanism/geo_ind`). Points can be added with ``+``.

   .. py:attribute:: x
      :type: float

   .. py:attribute:: y
      :type: float

      The coordinates of the point.

   .. py:staticmethod:: from_polar(radius, angle)

      The point with the given polar coordinates (``angle`` in radians).

   .. py:staticmethod:: from_cell(cell, width)

      The point ``(cell % width, cell // width)``: the coordinates of the cell with index ``cell`` in a grid with
      ``width`` columns (see :py:func:`qif.probab.to_grid`).
```
