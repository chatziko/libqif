import math

import pytest
from qif.measure import bayes_vuln

from helpers import assert_equal, num


def test_vulnerability(eT, t):
	assert_equal(num(eT, "5/10"), bayes_vuln.prior(t.unif_2))
	assert_equal(num(eT, "1/10"), bayes_vuln.prior(t.unif_10))
	assert_equal(num(eT, 1),      bayes_vuln.prior(t.point_4))
	assert_equal(num(eT, "8/10"), bayes_vuln.prior(t.pi1))

def test_post_vulnerability(eT, t):
	assert_equal(num(eT, 1),      bayes_vuln.posterior(t.unif_2, t.id_2))
	assert_equal(num(eT, 1),      bayes_vuln.posterior(t.point_2, t.id_2))
	assert_equal(num(eT, 1),      bayes_vuln.posterior(t.pi1, t.id_2))

	assert_equal(num(eT, 1),      bayes_vuln.posterior(t.unif_10, t.id_10))
	assert_equal(num(eT, 1),      bayes_vuln.posterior(t.point_10, t.id_10))
	assert_equal(num(eT, 1),      bayes_vuln.posterior(t.pi2, t.id_10))

	assert_equal(num(eT, "1/10"), bayes_vuln.posterior(t.unif_10, t.noint_10))
	assert_equal(num(eT, 1),      bayes_vuln.posterior(t.point_10, t.noint_10))

	assert_equal(bayes_vuln.prior(t.pi2), bayes_vuln.posterior(t.pi2, t.noint_10))

	assert_equal(bayes_vuln.prior(t.pi3), bayes_vuln.posterior(t.pi3, t.c1))	# no change in vulnerability
	assert_equal(num(eT, "31/40"), bayes_vuln.posterior(t.pi4, t.c1))

	with pytest.raises(RuntimeError, match="invalid prior size"):
		bayes_vuln.posterior(t.unif_2, t.id_10)

def test_mult_capacity(eT, t):
	assert_equal(num(eT, 2),       bayes_vuln.mult_capacity(t.id_2))
	assert_equal(num(eT, 10),      bayes_vuln.mult_capacity(t.id_10))
	assert_equal(num(eT, 1),       bayes_vuln.mult_capacity(t.noint_10))
	assert_equal(num(eT, "15/10"), bayes_vuln.mult_capacity(t.c1))

	# max_mult_leakage is given for uniform prior
	assert_equal(bayes_vuln.mult_leakage(t.unif_10, t.crand_10), bayes_vuln.mult_capacity(t.crand_10))


def test_min_entropy_leakage(real_eT, rt):
	assert_equal(1,              bayes_vuln.min_entropy_leakage(rt.unif_2, rt.id_2))
	assert_equal(math.log2(10),  bayes_vuln.min_entropy_leakage(rt.unif_10, rt.id_10))
	assert_equal(0,              bayes_vuln.min_entropy_leakage(rt.unif_10, rt.noint_10))
	assert_equal(math.log2(1.5), bayes_vuln.min_entropy_leakage(rt.unif_2, rt.c1))

def test_mult_capacity_bound_cap(real_eT, rt):
	assert_equal(627991708.193414211273193359375, bayes_vuln.mult_capacity_bound_cap(rt.id_4, 1000000))
