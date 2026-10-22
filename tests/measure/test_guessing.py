import pytest
from qif.measure import guessing

from helpers import assert_equal, num


def test_vulnerability(eT, t):
	assert_equal(num(eT, "15/10"), guessing.prior(t.unif_2))
	assert_equal(num(eT, "55/10"), guessing.prior(t.unif_10))
	assert_equal(num(eT, 1),       guessing.prior(t.point_4))
	assert_equal(num(eT, "12/10"), guessing.prior(t.pi1))

def test_post_entropy(eT, t):
	assert_equal(num(eT, 1),       guessing.posterior(t.unif_2, t.id_2))
	assert_equal(num(eT, 1),       guessing.posterior(t.point_2, t.id_2))
	assert_equal(num(eT, 1),       guessing.posterior(t.pi1, t.id_2))

	assert_equal(num(eT, 1),       guessing.posterior(t.unif_10, t.id_10))
	assert_equal(num(eT, 1),       guessing.posterior(t.point_10, t.id_10))
	assert_equal(num(eT, 1),       guessing.posterior(t.pi2, t.id_10))

	assert_equal(num(eT, "55/10"), guessing.posterior(t.unif_10, t.noint_10))
	assert_equal(num(eT, 1),       guessing.posterior(t.point_10, t.noint_10))

	assert_equal(guessing.prior(t.pi2), guessing.posterior(t.pi2, t.noint_10))

	assert_equal(guessing.prior(t.pi3), guessing.posterior(t.pi3, t.c1))	# no change in entropy
	assert_equal(num(eT, "1225/1000"), guessing.posterior(t.pi4, t.c1))

	with pytest.raises(RuntimeError):
		guessing.posterior(t.unif_2, t.id_10)
