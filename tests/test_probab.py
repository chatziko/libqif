import pytest
from qif import probab

from helpers import arr, assert_arr_equal, assert_proper_prob


def test_construct(eT):
	pi = arr(eT, "0.5 0.25 0.25")
	probab.assert_proper(pi)

	# malformed prob
	with pytest.raises(RuntimeError, match="not a proper dist"):
		probab.assert_proper(arr(eT, "0.1 0.2 0.3"))

def test_uniform(eT):
	assert_proper_prob(probab.uniform(1), 1)
	assert_arr_equal(probab.uniform(4), arr(eT, "0.25 0.25 0.25 0.25"))

def test_randu(eT):
	assert_proper_prob(probab.randu(200), 200)
	assert_proper_prob(probab.randu(5), 5)

def test_point(eT):
	assert_arr_equal(probab.point(4),    arr(eT, "1 0 0 0"))
	assert_arr_equal(probab.point(4, 2), arr(eT, "0 0 1 0"))
