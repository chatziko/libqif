"""
Quantitative Information Flow library.

.. autosummary::
	:toctree: _autosummary
	:template: template.rst

	channel
	probab
	metric
	measure
	mechanism
	refinement
	utility
	lp

|
"""

try:
	import numpy
except:
	raise Exception("numpy is required by qif")

# In the cpp code, "_qif" is not a package, although it contains other sub-modules! Same
# for qif.channel, qif.metric, qif.measure. Note that "import _qif.channel" works, although
# _qif is not a package, simpy because '_qif.channel' is stored in sys.modules when _qif loads.
# 
# To make the code cleaner, we follow the standard python approach and create packages for
# these, and import their contents in the corresponding __init__.py.
#
from . import channel, metric, measure, mechanism			# packages
from ._qif import probab, refinement, utility, lp			# modules
from ._qif import __version__, point, set_default_type		# other stuff

# data type aliases
from numpy import float64 as double
from fractions import Fraction as rat
from numpy import uint32 as uint

# numpy formatter, so that rats are nicely displayed
numpy.set_printoptions(formatter={ "object": str })
