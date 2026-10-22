from qif import metric
from qif.measure import d_privacy as measure_dp
from qif.mechanism import geo_ind

from helpers import assert_equal, assert_proper_chan, grid


def test_grid(real_eT):
	width = height = 3
	size = width * height
	step = 1.1
	epsilon = 0.9

	d = grid(width, step)

	laplace = geo_ind.planar_laplace_grid(width, height, step, epsilon)
	geom = geo_ind.planar_geometric_grid(width, height, step, epsilon)

	# for the planar laplace, the accuracy that we get through numeric integration is not that great
	assert_proper_chan(laplace, size, size, mrd=1e-3)
	assert measure_dp.is_private(laplace, metric.scale(d, epsilon))
	assert_equal(0.8844, measure_dp.smallest_epsilon(laplace, d), md=0, mrd=1e-3)

	# same (just a bit better) for planar_geometric
	assert_proper_chan(geom, size, size)
	assert measure_dp.is_private(geom, metric.scale(d, epsilon + 1e-5))
	assert_equal(epsilon, measure_dp.smallest_epsilon(geom, d), md=0, mrd=1e-6)
