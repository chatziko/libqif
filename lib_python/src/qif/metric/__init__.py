"""
Metrics and other distance functions.
"""

import builtins
import math

from .._qif.metric import *

# The metric combinators are implemented in python, so that they work for metrics on any domain (secrets, numbers,
# points, distributions). C++ bindings would need to fix the domain, since pybind11 wraps any python callable passed
# as std::function<R(T,T)>, converting the arguments to T.
#
# No docstrings: the documentation is in misc/docs.

def scale(d, coeff):
	# 0 stays 0 even if coeff is inf
	return lambda a, b: (lambda r: r * coeff if r != 0 else r)(d(a, b))

def min(d1, d2):
	return lambda a, b: builtins.min(d1(a, b), d2(a, b))

def max(d1, d2):
	return lambda a, b: builtins.max(d1(a, b), d2(a, b))

def mirror(d):
	return lambda a, b: d(b, a)

def threshold(d, thres):
	return lambda a, b: builtins.min(d(a, b), thres)

def threshold_inf(d, thres):
	return lambda a, b: (lambda r: r if r <= thres else math.inf)(d(a, b))
