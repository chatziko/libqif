# Helpers for the qif tests (ported from tests_cpp/include/tests_aux.h)

from fractions import Fraction
from functools import cached_property

import numpy
from qif import channel, double, metric, point, probab, rat


def num(eT, x):
	"""Number of type eT, from an int, a "0.8" or "4/10" string, or a Fraction."""
	return Fraction(x) if eT is rat else double(Fraction(x))

def arr(eT, s):
	"""Vector/matrix of type eT from a matlab-like string, eg "0.8 0.2; 0.3 0.7"."""
	rows = [[num(eT, x) for x in row.split()] for row in s.split(";") if row.strip()]
	return numpy.array(rows[0] if len(rows) == 1 else rows, dtype=object if eT is rat else double)


# Same semantics as qif::equal in lib_cpp/include/qif_bits/lib_aux.h: exact for rat; for double, absolute
# error md for comparisons with 0, otherwise mrd relative to the largest number (as absolute error up to 1).
DEF_MD = 1e-7
DEF_MRD = 100 * numpy.finfo(double).eps

def equal(x, y, md=DEF_MD, mrd=DEF_MRD):
	if isinstance(x, Fraction) or isinstance(y, Fraction):
		# both should be rat, so that a computation silently done in double is detected (0.5 == Fraction(1,2) in python)
		return isinstance(x, Fraction) and isinstance(y, Fraction) and x == y
	x, y = float(x), float(y)
	if x == y:
		return True
	largest = max(abs(x), abs(y))
	if largest == float("inf"):
		return False
	return abs(x - y) <= (md if x == 0 or y == 0 or mrd == 0 else mrd * max(largest, 1))

def assert_equal(x, y, **kw):
	assert equal(x, y, **kw), f"{x} and {y} are not equal"

def assert_arr_equal(A, B, **kw):
	"""Element-wise equal (same shape), for channels and distributions."""
	A, B = numpy.asarray(A), numpy.asarray(B)
	assert A.shape == B.shape and all(equal(x, y, **kw) for x, y in zip(A.flat, B.flat)), f"not equal:\n{A}\n{B}"

def assert_proper_chan(C, n_rows, n_cols, **kw):
	assert C.shape == (n_rows, n_cols), f"shape {C.shape} != {(n_rows, n_cols)}"
	assert channel.is_proper(C, **kw), f"not a proper channel:\n{C}"

def assert_proper_prob(pi, n):
	assert pi.shape == (n,), f"shape {pi.shape} != {(n,)}"
	assert probab.is_proper(pi), f"not a proper distribution: {pi}"


def grid(width, step=1):
	"""Euclidean metric on the cells of a grid of the given width (C++'s metric::grid, scaled by step)."""
	euclid = metric.euclidean(point)
	return lambda a, b: step * euclid(point.from_cell(a, width), point.from_cell(b, width))


class Base:
	"""Commonly used distributions/channels of the default type (BaseTest in tests_aux.h). Created on first access."""

	def __init__(self, eT):
		self.eT = eT

	unif_2   = cached_property(lambda t: probab.uniform(2))
	unif_4   = cached_property(lambda t: probab.uniform(4))
	unif_10  = cached_property(lambda t: probab.uniform(10))
	unif_100 = cached_property(lambda t: probab.uniform(100))
	point_2  = cached_property(lambda t: probab.point(2))
	point_4  = cached_property(lambda t: probab.point(4))
	point_10 = cached_property(lambda t: probab.point(10))
	pi1      = cached_property(lambda t: arr(t.eT, "0.2 0.8"))
	pi2      = cached_property(lambda t: arr(t.eT, "0.2 0.8 0 0 0 0 0 0 0 0"))
	pi3      = cached_property(lambda t: arr(t.eT, "0.25 0.75"))
	pi4      = cached_property(lambda t: arr(t.eT, "0.75 0.25"))
	pi5      = cached_property(lambda t: arr(t.eT, "0.1 0.1 0.1 0.7"))
	prand_10 = cached_property(lambda t: probab.randu(10))

	id_2      = cached_property(lambda t: channel.identity(2))
	id_4      = cached_property(lambda t: channel.identity(4))
	id_10     = cached_property(lambda t: channel.identity(10))
	noint_4   = cached_property(lambda t: channel.no_interference(4))
	noint_10  = cached_property(lambda t: channel.no_interference(10))
	c1        = cached_property(lambda t: arr(t.eT, "0.8 0.2; 0.3 0.7"))
	crand_2   = cached_property(lambda t: channel.randu(2))
	crand_10  = cached_property(lambda t: channel.randu(10))
	crand_100 = cached_property(lambda t: channel.randu(100))
