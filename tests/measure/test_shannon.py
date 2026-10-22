import math

import pytest
from qif import probab
from qif.measure import shannon

from helpers import arr, assert_arr_equal, assert_equal


def test_entropy(real_eT, rt):
	assert_equal(1,                 shannon.prior(rt.unif_2))
	assert_equal(math.log2(10),     shannon.prior(rt.unif_10))
	assert_equal(0,                 shannon.prior(rt.point_4))
	assert_equal(0.721928094887362, shannon.prior(rt.pi1))

def test_cond_entropy(real_eT, rt):
	assert_equal(0, shannon.posterior(rt.unif_2, rt.id_2))
	assert_equal(0, shannon.posterior(rt.point_2, rt.id_2))
	assert_equal(0, shannon.posterior(rt.pi1, rt.id_2))

	assert_equal(0, shannon.posterior(rt.unif_10, rt.id_10))
	assert_equal(0, shannon.posterior(rt.point_10, rt.id_10))
	assert_equal(0, shannon.posterior(rt.pi2, rt.id_10))

	assert_equal(math.log2(10), shannon.posterior(rt.unif_10, rt.noint_10))
	assert_equal(0,             shannon.posterior(rt.point_10, rt.noint_10))

	pi = probab.randu(10)
	assert_equal(shannon.prior(pi), shannon.posterior(pi, rt.noint_10))

	assert_equal(0.669020059980807, shannon.posterior(rt.pi3, rt.c1))

	with pytest.raises(RuntimeError):
		shannon.posterior(rt.unif_2, rt.id_10)

def test_capacity(real_eT, rt):
	cap1, pi1 = shannon.add_capacity(rt.id_2)
	assert_equal(1, cap1)
	assert_arr_equal(pi1, rt.unif_2)

	cap2, pi2 = shannon.add_capacity(rt.id_10)
	assert_equal(math.log2(10), cap2)
	assert_arr_equal(pi2, rt.unif_10)

	assert_equal(0, shannon.add_capacity(rt.noint_10)[0], mrd=0)

	assert_equal(0.19123813831431799, shannon.add_capacity(rt.c1)[0])

	# symmetric
	C = arr(real_eT,
		".3 .2 .5;"
		".5 .3 .2;"
		".2 .5 .3;"
	)
	cap = math.log2(C.shape[1]) - shannon.prior(C[0])
	assert_equal(cap, shannon.add_capacity(C)[0])

	# weakly symmetric
	C = arr(real_eT,
		"0.333333333 0.166666667 0.5;"
		"0.333333333 0.5         0.166666667;"
	)
	cap = math.log2(C.shape[1]) - shannon.prior(C[0])
	assert_equal(cap, shannon.add_capacity(C)[0])
