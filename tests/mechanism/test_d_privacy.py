import math

import pytest
import qif
from qif import channel, metric, probab, uint
from qif.measure import d_privacy as measure_dp
from qif.mechanism import d_privacy

from helpers import assert_arr_equal, assert_equal, assert_proper_chan, grid


def test_reals(real_eT):
	size = 10
	step = 1.1
	epsilon = 0.9

	d = metric.scale(metric.euclidean(uint), step)

	geom = d_privacy.geometric(size, epsilon * step)
	expon = d_privacy.exponential(size, metric.scale(d, epsilon))
	tc = d_privacy.tight_constraints(size, metric.scale(d, epsilon))

	assert_proper_chan(geom, size, size)
	assert     measure_dp.is_private(geom, metric.scale(d, epsilon))
	assert not measure_dp.is_private(geom, metric.scale(d, epsilon - 0.01))
	assert_equal(epsilon, measure_dp.smallest_epsilon(geom, d))

	assert_proper_chan(expon, size, size)
	assert measure_dp.is_private(expon, metric.scale(d, epsilon))
	assert_equal(0.64197307180467134, measure_dp.smallest_epsilon(expon, d))

	assert_arr_equal(geom, tc)

	# fat and skinny
	for n_rows, n_cols in [(size, 2*size), (2*size, size)]:
		geom = d_privacy.geometric(n_rows, epsilon * step, n_cols)
		expon = d_privacy.exponential(n_rows, metric.scale(d, epsilon), n_cols)

		assert_proper_chan(geom, n_rows, n_cols)
		assert     measure_dp.is_private(geom, metric.scale(d, epsilon))
		assert not measure_dp.is_private(geom, metric.scale(d, epsilon - 0.01))
		assert_equal(epsilon, measure_dp.smallest_epsilon(geom, d))

		assert_proper_chan(expon, n_rows, n_cols)
		assert measure_dp.is_private(expon, metric.scale(d, epsilon))

def test_geometric_truncation(real_eT):
	# We create 4 channels with secrets 5..14 and obs 0..9 / 10..19 / 1..18 / 7..12
	# We first create by remapping a 0..19 x 0..19 channel, then create directly and compare
	#
	epsilon = 0.9
	G0 = d_privacy.geometric(20, epsilon)		# the 0..19 x 0..19 to remap

	clamp = lambda lo, hi: channel.deterministic(lambda x: min(max(x, lo), hi), 20, 20)
	T1 = G0 @ clamp(0,  9)			# 0..9
	T2 = G0 @ clamp(10, 19)			# 10..19
	T3 = G0 @ clamp(1,  18)			# 1..18
	T4 = G0 @ clamp(7,  12)			# 7..12

	# crop rows to 5..14 and the columns to 0..9 / 10..19 / 1..18 / 7..12,
	# compare to the same channels created directly by passing first_x / first_y
	assert_arr_equal(T1[5:15, 0:10],  d_privacy.geometric(10, epsilon, 10, 5, 0))
	assert_arr_equal(T2[5:15, 10:20], d_privacy.geometric(10, epsilon, 10, 5, 10))
	assert_arr_equal(T3[5:15, 1:19],  d_privacy.geometric(10, epsilon, 18, 5, 1))
	assert_arr_equal(T4[5:15, 7:13],  d_privacy.geometric(10, epsilon, 6,  5, 7))

def test_discrete(real_eT):
	size = 5
	epsilon = 0.9

	d = metric.discrete(uint)

	tc = d_privacy.tight_constraints(size, metric.scale(d, epsilon))
	expon = d_privacy.exponential(size, metric.scale(d, 2 * epsilon))
	rr = d_privacy.randomized_response(size, epsilon)

	assert_proper_chan(tc, size, size)
	assert     measure_dp.is_private(tc, metric.scale(d, epsilon))
	assert not measure_dp.is_private(tc, metric.scale(d, epsilon - 0.01))
	assert_equal(epsilon, measure_dp.smallest_epsilon(tc, d))

	# the exponential mechanism with 2*epsilon should be the same as the tight constraints and randomized response
	assert_arr_equal(expon, tc)
	assert_arr_equal(expon, rr)

def test_grid(real_eT):
	width = height = 3
	size = width * height
	step = 1.1
	epsilon = 0.9

	d = grid(width, step)

	tc = d_privacy.tight_constraints(size, metric.scale(d, epsilon))
	expon = d_privacy.exponential(size, metric.scale(d, epsilon))

	assert_proper_chan(tc, size, size)
	assert     measure_dp.is_private(tc, metric.scale(d, epsilon))
	assert not measure_dp.is_private(tc, metric.scale(d, epsilon - 0.01))
	assert_equal(epsilon, measure_dp.smallest_epsilon(tc, d))

	assert_proper_chan(expon, size, size)
	assert measure_dp.is_private(expon, metric.scale(d, epsilon))
	assert_equal(0.58945591528726249, measure_dp.smallest_epsilon(expon, d))

def test_opt_exp_loss(real_eT, rt):
	size = 10
	step = 1.1
	epsilon = 0.9

	# with euclidean metric we get the geometric mechanism
	d = metric.scale(metric.euclidean(uint), step)
	bin_loss = metric.discrete(uint)

	geom = d_privacy.geometric(size, epsilon * step)
	opt = d_privacy.min_loss_given_d(rt.unif_10, size, metric.scale(d, epsilon), bin_loss)
	dist_opt = d_privacy.min_loss_given_d(rt.unif_10, size, metric.scale(d, epsilon), bin_loss, "dist")
	dist_opt2 = d_privacy.min_loss_given_d(rt.unif_10, size, metric.scale(d, epsilon), bin_loss, "dist_strict")

	assert_arr_equal(geom, opt)
	assert_arr_equal(geom, dist_opt, md=1e-5, mrd=0)		# we get good accuracy on linux but not OSX!

	assert channel.is_proper(dist_opt2, 1e-3)				# not sure why we don't get good precision

	# with inf constraints we get the identity channel
	d = metric.scale(metric.discrete(uint), math.inf)
	opt = d_privacy.min_loss_given_d(rt.unif_10, size, metric.scale(d, epsilon), bin_loss)

	assert_arr_equal(rt.id_10, opt)

def test_opt_exp_loss_badly_scaled(real_eT):
	# large epsilon gives coefficients up to exp(380), CLP/HiGHS/GLOP fail on this program
	size = 20
	d = metric.scale(metric.euclidean(uint), 20)
	loss = metric.euclidean(uint)

	# AUTO retries with SoPlex when the chosen solver fails
	opt = d_privacy.min_loss_given_d(probab.uniform(size), size, d, loss)
	assert_proper_chan(opt, size, size)

	# an explicitly chosen solver that fails gives an error (instead of an empty channel)
	qif.lp.defaults.solver = "CLP"
	try:
		with pytest.raises(RuntimeError, match="failed|not available"):	# not available: built without OR-Tools
			d_privacy.min_loss_given_d(probab.uniform(size), size, d, loss)
	finally:
		qif.lp.defaults.solver = "AUTO"
