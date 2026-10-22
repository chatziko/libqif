import pytest
from qif.measure import g_vuln

from helpers import assert_equal, num

# TODO: test more gain functions. Currently these tests are copies from test_bayes_vuln.py and test only the id gain function


def test_vulnerability(eT, t):
	assert_equal(num(eT, "5/10"), g_vuln.prior(t.id_2, t.unif_2))
	assert_equal(num(eT, "1/10"), g_vuln.prior(t.id_10, t.unif_10))
	assert_equal(num(eT, 1),      g_vuln.prior(t.id_4, t.point_4))
	assert_equal(num(eT, "8/10"), g_vuln.prior(t.id_2, t.pi1))

	# the same with the g_id gain function (instead of the G matrix)
	assert_equal(num(eT, "1/10"), g_vuln.prior(g_vuln.g_id(), t.unif_10))
	assert_equal(num(eT, "8/10"), g_vuln.prior(g_vuln.g_id(), t.pi1))

def test_post_vulnerability(eT, t):
	assert_equal(num(eT, 1),      g_vuln.posterior(t.id_2, t.unif_2, t.id_2))
	assert_equal(num(eT, 1),      g_vuln.posterior(t.id_2, t.point_2, t.id_2))
	assert_equal(num(eT, 1),      g_vuln.posterior(t.id_2, t.pi1, t.id_2))

	assert_equal(num(eT, 1),      g_vuln.posterior(t.id_10, t.unif_10, t.id_10))
	assert_equal(num(eT, 1),      g_vuln.posterior(t.id_10, t.point_10, t.id_10))
	assert_equal(num(eT, 1),      g_vuln.posterior(t.id_10, t.pi2, t.id_10))

	assert_equal(num(eT, "1/10"), g_vuln.posterior(t.id_10, t.unif_10, t.noint_10))
	assert_equal(num(eT, 1),      g_vuln.posterior(t.id_10, t.point_10, t.noint_10))

	assert_equal(g_vuln.prior(t.id_10, t.pi2), g_vuln.posterior(t.id_10, t.pi2, t.noint_10))

	assert_equal(g_vuln.prior(t.id_2, t.pi3), g_vuln.posterior(t.id_2, t.pi3, t.c1))	# no change in entropy
	assert_equal(num(eT, "31/40"), g_vuln.posterior(t.id_2, t.pi4, t.c1))
	assert_equal(num(eT, "31/40"), g_vuln.posterior(g_vuln.g_id(), t.pi4, t.c1))

	with pytest.raises(RuntimeError):
		g_vuln.posterior(t.id_10, t.unif_2, t.id_10)

def test_add_capacity(eT, t):
	assert_equal(num(eT, "1/2"), g_vuln.add_capacity(t.unif_2,  t.id_2,    True))
	assert_equal(num(eT, 1),     g_vuln.add_capacity(t.unif_2,  t.id_2,    False))

	assert_equal(num(eT, 0),     g_vuln.add_capacity(t.point_2, t.id_2,    True))
	assert_equal(num(eT, 0),     g_vuln.add_capacity(t.point_2, t.id_2,    False))

	assert_equal(num(eT, 0),     g_vuln.add_capacity(t.unif_4,  t.noint_4, True))
	assert_equal(num(eT, 0),     g_vuln.add_capacity(t.unif_4,  t.noint_4, False))

	with pytest.raises(RuntimeError):
		g_vuln.add_capacity(t.unif_2, t.id_10)
