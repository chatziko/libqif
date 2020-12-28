import sys
from typing import (
	Callable as Callable,
	Tuple as Tuple,
	overload as overload,
	Union,
	Type,
	TypeVar,
	Any as Any,
	List as List,
)
from numpy import ndarray as ndarray, array as array, float64 as double, uint32 as uint
from fractions import Fraction as rat
from . import point as point

TypeLike = Union[Type[double], Type[rat], Type[point], Type[uint]]
FloatOrRat = Union[float, rat]

R = TypeVar('R')
T = TypeVar('T')

Metric = Callable[[T, T], R]

def_type: TypeLike
def_md: FloatOrRat
def_mrd: FloatOrRat
inf: FloatOrRat
