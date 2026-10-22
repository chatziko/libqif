import math

import pytest
from qif import metric, point, probab, rat, uint

from helpers import DEF_MRD, arr, assert_arr_equal, assert_equal, num

# Note: metrics needing log/sqrt are only available for double, the rest are tested for both types


def test_euclidean_uint(eT):
	# distances on secrets are of the default type
	euclid = metric.euclidean(uint)
	chainable = metric.euclidean_chain(uint)

	assert_equal(num(eT, 0), euclid(3, 3))
	assert_equal(num(eT, 5), euclid(0, 5))
	assert_equal(num(eT, 5), euclid(5, 0))

	assert not chainable(0, 1)
	assert not chainable(5, 6)
	assert     chainable(0, 5)

def test_euclidean_rat(eT):
	# for rat (the default type) the euclidean metric is on rat values
	euclid = metric.euclidean()

	assert_equal(eT(0), euclid(eT(3), eT(3)))
	assert_equal(eT(5), euclid(eT(0), eT(5)))
	assert_equal(eT(5), euclid(eT(5), eT(0)))

def test_discrete(eT):
	disc = metric.discrete(uint)

	assert_equal(num(eT, 0), disc(3, 3))
	assert_equal(num(eT, 1), disc(0, 2))
	assert_equal(num(eT, 1), disc(5, 0))

def test_scale(eT):
	scaled_euclid = metric.scale(metric.euclidean(uint), num(eT, 10))
	scaled_disc   = metric.scale(metric.discrete(uint), num(eT, 10))

	assert_equal(num(eT, 0),  scaled_euclid(3, 3))
	assert_equal(num(eT, 50), scaled_euclid(0, 5))
	assert_equal(num(eT, 50), scaled_euclid(5, 0))

	assert_equal(num(eT, 0),  scaled_disc(3, 3))
	assert_equal(num(eT, 10), scaled_disc(0, 3))
	assert_equal(num(eT, 10), scaled_disc(5, 0))

def test_scale_inf(real_eT):
	# 0 stays 0 even if coeff is inf (rat has no inf)
	assert_equal(0,        metric.scale(metric.discrete(uint), math.inf)(3, 3))
	assert_equal(math.inf, metric.scale(metric.discrete(uint), math.inf)(0, 3))

def test_threshold(real_eT):
	thres = metric.threshold(metric.euclidean(uint), 10)
	thres_inf = metric.threshold_inf(metric.euclidean(uint), 10)

	assert_equal(0,  thres(3, 3))
	assert_equal(5,  thres(0, 5))
	assert_equal(10, thres(0, 10))
	assert_equal(10, thres(0, 20))

	assert_equal(5,        thres_inf(0, 5))
	assert_equal(10,       thres_inf(0, 10))
	assert_equal(math.inf, thres_inf(0, 20))

def test_euclidean_point(real_eT):
	euclid = metric.euclidean(point)

	assert_equal(0,            euclid(point(1, 1), point(1, 1)))
	assert_equal(math.sqrt(2), euclid(point(0, 0), point(1, 1)))
	assert_equal(math.sqrt(5), euclid(point(2, 3), point(1, 1)))

def test_manhattan_point(real_eT):
	manh = metric.manhattan()

	assert_equal(0, manh(point(1, 1), point(1, 1)))
	assert_equal(2, manh(point(0, 0), point(1, 1)))
	assert_equal(3, manh(point(2, 3), point(1, 1)))

def test_grid_point(real_eT):
	# euclidean/manhattan on the cells of a 4-wide grid
	euclid, manh = metric.euclidean(point), metric.manhattan()
	grid_euclid = lambda a, b: euclid(point.from_cell(a, 4), point.from_cell(b, 4))
	grid_manh   = lambda a, b: manh  (point.from_cell(a, 4), point.from_cell(b, 4))

	assert_equal(0,            grid_euclid(5, 5))
	assert_equal(math.sqrt(2), grid_euclid(0, 5))		# cell  5 is (1,1)
	assert_equal(math.sqrt(5), grid_euclid(11, 5))		# cell 11 is (3,2)

	assert_equal(0, grid_manh(5, 5))
	assert_equal(2, grid_manh(0, 5))
	assert_equal(3, grid_manh(11, 5))

def test_l1(eT, t):
	l1 = metric.l1()

	assert_equal(num(eT, 0),     l1(t.unif_4, t.unif_4))
	assert_equal(num(eT, "3/2"), l1(t.unif_4, t.point_4))
	assert_equal(num(eT, "11/10"), l1(t.pi1, arr(eT, "0.75 0.25")))	# |0.2-0.75| + |0.8-0.25|

	# rat should give exact results, not fall back to double (0.75 == Fraction(3,4) so check the type)
	if eT is rat:
		assert type(l1(t.unif_4, t.point_4)) is rat

def test_total_variation(eT, t):
	tv = metric.total_variation()

	assert_equal(num(eT, 0),     tv(t.unif_4, t.unif_4))
	assert_equal(num(eT, 0),     tv(t.point_4, t.point_4))
	assert_equal(num(eT, "3/4"), tv(t.unif_4, t.point_4))
	assert_equal(num(eT, "3/4"), tv(t.point_4, t.unif_4))

	if eT is rat:
		assert type(tv(t.unif_4, t.point_4)) is rat

def test_convex_separation(eT, t):
	bed = metric.convex_separation()

	assert_equal(num(eT, 0),      bed(t.unif_4, t.unif_4))
	assert_equal(num(eT, 0),      bed(t.point_4, t.point_4))

	assert_equal(num(eT, 1),      bed(t.unif_4, t.point_4))
	assert_equal(num(eT, 1),      bed(t.point_4, t.unif_4))

	assert_equal(num(eT, "9/14"), bed(t.unif_4, t.pi5))
	assert_equal(num(eT, "9/14"), bed(t.pi5, t.unif_4))

	# quasi: not symmetric
	q = metric.convex_separation_quasi()
	assert_equal(num(eT, 1),      q(t.unif_4, t.point_4))
	assert_equal(num(eT, "3/4"),  q(t.point_4, t.unif_4))

def test_distance_matrix(eT):
	M = metric.to_distance_matrix(metric.euclidean(uint), 3)
	assert_arr_equal(M, arr(eT, "0 1 2; 1 0 1; 2 1 0"))

	d = metric.from_distance_matrix(M)
	assert_equal(num(eT, 2), d(0, 2))
	assert_equal(num(eT, 1), d(2, 1))

def test_multiplicative_distance(real_eT, rt):
	mtv = metric.mult_total_variation()

	assert_equal(0,                  mtv(rt.unif_4, rt.unif_4))
	assert_equal(0,                  mtv(rt.point_4, rt.point_4))

	assert_equal(math.inf,           mtv(rt.unif_4, rt.point_4))
	assert_equal(math.inf,           mtv(rt.point_4, rt.unif_4))

	assert_equal(math.log(0.7/0.25), mtv(rt.unif_4, rt.pi5))
	assert_equal(math.log(0.7/0.25), mtv(rt.pi5, rt.unif_4))

@pytest.mark.parametrize("method", ["auto", "fastemd", "lp"])
def test_kantorovich(eT, t, method):
	if method == "fastemd" and eT is rat:
		with pytest.raises(RuntimeError, match="fastemd method only available for double"):
			metric.kantorovich(metric.discrete(uint), method)
		return

	disc        = metric.discrete(uint)
	euclid      = metric.euclidean(uint)
	kant_disc   = metric.kantorovich(disc, method)
	kant_euclid = metric.kantorovich(euclid, method)
	tv          = metric.total_variation()

	# fastemd ("auto" for double) has worse precision, so we need a larger mrd
	mrd = 1e-5 if method != "lp" and eT is not rat else DEF_MRD

	assert_equal(num(eT, 0),     kant_disc(t.unif_4, t.unif_4))
	assert_equal(num(eT, 0),     kant_disc(t.point_4, t.point_4))
	assert_equal(num(eT, "3/4"), kant_disc(t.unif_4, t.point_4))
	assert_equal(num(eT, "3/4"), kant_disc(t.point_4, t.unif_4))

	assert_equal(num(eT, 0),     kant_euclid(t.unif_4, t.unif_4))
	assert_equal(num(eT, 0),     kant_euclid(t.point_4, t.point_4))
	assert_equal(num(eT, "3/2"), kant_euclid(t.unif_4, t.point_4))
	assert_equal(num(eT, "3/2"), kant_euclid(t.point_4, t.unif_4))
	assert_equal(num(eT, "9/2"), kant_euclid(t.unif_10, t.point_10))
	assert_equal(num(eT, "9/2"), kant_euclid(t.point_10, t.unif_10))

	for i in range(3):
		# distance of point dists is the same as the distance between the corresponding elements
		assert_equal(disc  (0, i), kant_disc  (t.point_4, probab.point(4, i)))
		assert_equal(euclid(0, i), kant_euclid(t.point_4, probab.point(4, i)), md=0, mrd=mrd)

	# kantorovich over the discrete metric = total variation
	# (only for double: for rat the random distributions have huge denominators, too slow)
	if eT is not rat:
		for _ in range(10):
			p1, p2 = probab.randu(10), probab.randu(10)
			assert_equal(tv(p1, p2), kant_disc(p1, p2), md=0, mrd=mrd)

def test_kantorovich_fastemd_vs_lp(real_eT):
	# both methods should produce the same result
	disc, euclid = metric.discrete(uint), metric.euclidean(uint)

	for _ in range(10):
		p1, p2 = probab.randu(10), probab.randu(10)
		assert_equal(metric.kantorovich(disc, "fastemd")(p1, p2),   metric.kantorovich(disc, "lp")(p1, p2),   md=0, mrd=1e-5)
		assert_equal(metric.kantorovich(euclid, "fastemd")(p1, p2), metric.kantorovich(euclid, "lp")(p1, p2), md=0, mrd=1e-4)

	with pytest.raises(RuntimeError, match="invalid method: foo"):
		metric.kantorovich(disc, "foo")

def test_mult_kantorovich(real_eT, rt):
	disc         = metric.scale(metric.discrete(uint), math.inf)
	euclid       = metric.euclidean(uint)
	mkant_disc   = metric.mult_kantorovich(disc)
	mkant_euclid = metric.mult_kantorovich(euclid)
	mtv          = metric.mult_total_variation()

	assert_equal(0,                  mkant_disc(rt.unif_4,  rt.unif_4))
	assert_equal(0,                  mkant_disc(rt.point_4, rt.point_4))
	assert_equal(math.inf,           mkant_disc(rt.unif_4,  rt.point_4))
	assert_equal(math.inf,           mkant_disc(rt.point_4, rt.unif_4))
	assert_equal(math.log(0.7/0.25), mkant_disc(rt.unif_4,  rt.pi5))
	assert_equal(math.log(0.7/0.25), mkant_disc(rt.pi5,     rt.unif_4))

	assert_equal(0,                  mkant_euclid(rt.unif_4,   rt.unif_4))
	assert_equal(0,                  mkant_euclid(rt.point_4,  rt.point_4))
	assert_equal(2.0538953374413045, mkant_euclid(rt.unif_4,   rt.point_4))
	assert_equal(2.0538953374413045, mkant_euclid(rt.point_4,  rt.unif_4))
	assert_equal(7.156044651432,     mkant_euclid(rt.unif_10,  rt.point_10), md=1e-12, mrd=0)
	assert_equal(7.156044651432,     mkant_euclid(rt.point_10, rt.unif_10),  md=1e-12, mrd=0)

	for i in range(3):
		# distance of point dists is the same as the distance between the corresponding elements
		assert_equal(disc  (0, i), mkant_disc  (rt.point_4, probab.point(4, i)))
		assert_equal(euclid(0, i), mkant_euclid(rt.point_4, probab.point(4, i)))

	# mult kantorovich over the discrete metric (with inf value) = mult total variation
	p1, p2 = probab.randu(10), probab.randu(10)
	assert_equal(mtv(p1, p2), mkant_disc(p1, p2))
