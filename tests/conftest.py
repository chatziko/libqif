# Shared pytest fixtures for the qif tests.
#
# Tests that take the 'eT' fixture run once per element type (double, rat), with that type set
# as the default, so probab.uniform(2), channel.identity(4), etc, produce elements of that type.
# Tests that only make sense for floating point take 'real_eT' instead. The 't' / 'rt' fixtures
# give the common distributions/channels of helpers.Base for 'eT' / 'real_eT' respectively.

import pytest
import qif
from qif import double, rat

from helpers import Base


ALL_TYPES = [double, rat]
REAL_TYPES = [double]

def _type_id(eT):
	return "rat" if eT is rat else "double"

def _with_default_type(eT):
	qif.set_default_type(eT)
	yield eT
	qif.set_default_type(double)

@pytest.fixture(params=ALL_TYPES, ids=_type_id)
def eT(request):
	yield from _with_default_type(request.param)

@pytest.fixture(params=REAL_TYPES, ids=_type_id)
def real_eT(request):
	yield from _with_default_type(request.param)


@pytest.fixture
def t(eT):
	return Base(eT)

@pytest.fixture
def rt(real_eT):
	return Base(real_eT)
