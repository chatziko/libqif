import math

import pytest
from qif import rat
from qif.measure import bayes_risk

from helpers import assert_equal, num


def test_vulnerability(eT, t):
	assert_equal(num(eT, "5/10"), bayes_risk.prior(t.unif_2))
	assert_equal(num(eT, "9/10"), bayes_risk.prior(t.unif_10))
	assert_equal(num(eT, 0),      bayes_risk.prior(t.point_4))
	assert_equal(num(eT, "2/10"), bayes_risk.prior(t.pi1))

def test_post_vulnerability(eT, t):
	assert_equal(num(eT, 0),      bayes_risk.posterior(t.unif_2, t.id_2))
	assert_equal(num(eT, 0),      bayes_risk.posterior(t.point_2, t.id_2))
	assert_equal(num(eT, 0),      bayes_risk.posterior(t.pi1, t.id_2))

	assert_equal(num(eT, 0),      bayes_risk.posterior(t.unif_10, t.id_10))
	assert_equal(num(eT, 0),      bayes_risk.posterior(t.point_10, t.id_10))
	assert_equal(num(eT, 0),      bayes_risk.posterior(t.pi2, t.id_10))

	assert_equal(num(eT, "9/10"), bayes_risk.posterior(t.unif_10, t.noint_10))
	assert_equal(num(eT, 0),      bayes_risk.posterior(t.point_10, t.noint_10))

	assert_equal(bayes_risk.prior(t.pi2), bayes_risk.posterior(t.pi2, t.noint_10))

	assert_equal(bayes_risk.prior(t.pi3), bayes_risk.posterior(t.pi3, t.c1))	# no change in vulnerability
	assert_equal(num(eT, "9/40"), bayes_risk.posterior(t.pi4, t.c1))

	with pytest.raises(RuntimeError, match="invalid prior size"):
		bayes_risk.posterior(t.unif_2, t.id_10)

def test_mult_capacity(eT, t):
	# the capacity of channels that can be broken (risk 0) is infinite. rat has no inf, it gives the max int64 instead.
	inf = num(eT, 2**63 - 1) if eT is rat else math.inf

	assert_equal(inf,            bayes_risk.mult_capacity(t.id_2)[0])
	assert_equal(inf,            bayes_risk.mult_capacity(t.id_10)[0])
	assert_equal(num(eT, 1),     bayes_risk.mult_capacity(t.noint_10)[0])
	assert_equal(num(eT, 2),     bayes_risk.mult_capacity(t.c1)[0])

	# max_mult_leakage is given for 2-secret uniform prior
	assert_equal(bayes_risk.mult_leakage(t.unif_2, t.crand_2), bayes_risk.mult_capacity(t.crand_2)[0])
