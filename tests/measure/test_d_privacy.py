import math

from qif import metric, uint
from qif.measure import d_privacy

from helpers import assert_equal


def test_is_private(real_eT, rt):
	d = metric.euclidean(uint)

	assert     d_privacy.is_private(rt.noint_10, metric.scale(d, 0))
	assert not d_privacy.is_private(rt.id_10,    metric.scale(d, 10000))

	assert     d_privacy.is_private(rt.c1, metric.scale(d, math.log(7.0/2)))
	assert not d_privacy.is_private(rt.c1, metric.scale(d, math.log(6.9/2)))

def test_smallest_epsilon(real_eT, rt):
	d = metric.euclidean(uint)

	assert_equal(0,               d_privacy.smallest_epsilon(rt.noint_10, d))
	assert_equal(math.inf,        d_privacy.smallest_epsilon(rt.id_10, d))
	assert_equal(math.log(7.0/2), d_privacy.smallest_epsilon(rt.c1, d))
