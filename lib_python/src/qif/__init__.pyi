"""
Quantitative Information Flow library.

Submodules
- channel
- probab
- metric
- measure
- mechanism
- refinement
- utility
- lp
"""
from . import typing as t
# from __future__ import annotations

# Note: "as" is needed for the modules to appear in pylance's auto-complete
#
from . import ( # packages
    channel as channel,
    metric as metric,
    measure as measure,
    mechanism as mechanism
)
from . import ( # modules
    probab as probab,
    refinement as refinement,
    utility as utility,
    lp as lp
)

from numpy import float64 as double
from fractions import Fraction as rat
from numpy import uint32 as uint

class point:
    def __add__(self, rhs: point) -> point: ...
    def __init__(self, x: float, y: float) -> None: ...
    def __repr__(self) -> str: ...

    @property
    def x(self) -> float: ...
    @x.setter
    def x(self, x: float) -> None: ...
    @property
    def y(self) -> float: ...
    @y.setter
    def y(self, y: float) -> None: ...

    @staticmethod
    def from_cell(cell: int, width: int) -> point: ...
    @staticmethod
    def from_polar(angle: float, radius: float) -> point: ...

def set_default_type(type: t.TypeLike) -> None: ...

__version__: str