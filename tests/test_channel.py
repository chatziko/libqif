import pytest
from qif import channel, double, mechanism, rat

from helpers import arr, assert_arr_equal, assert_proper_chan


def test_construct(eT):
	channel.assert_proper(arr(eT, "1 0 0; 0 1 0"))

	# malformed channel
	with pytest.raises(RuntimeError, match="not a proper matrix"):
		channel.assert_proper(arr(eT, "1 2; 3 0"))

def test_identity(eT):
	assert_proper_chan(channel.identity(0), 0, 0)
	assert_arr_equal(channel.identity(3), arr(eT, "1 0 0; 0 1 0; 0 0 1"))

def test_randu(eT):
	assert_proper_chan(channel.randu(200), 200, 200)
	assert_proper_chan(channel.randu(5), 5, 5)
	assert_proper_chan(channel.randu(4, 6), 4, 6)

@pytest.mark.parametrize("method", ["auto", "solve_lp", "project_subgrad", "project_lp", "project_qp", "lp", "subgrad"])	# lp, subgrad: deprecated names
def test_factorize(eT, t, method):
	if method in ("project_subgrad", "project_qp", "subgrad") and eT is rat:
		with pytest.raises(RuntimeError, match="not available for rat"):
			channel.factorize(t.id_4, t.noint_4, method=method)
		return

	# non factorizable
	assert_proper_chan(channel.factorize(t.id_10, t.noint_10, method=method), 0, 0)
	assert_proper_chan(channel.factorize(t.id_4,  t.noint_10, method=method), 0, 0)

	n, m = 4, 6
	B = channel.randu(n, m)
	A = B @ channel.randu(m, n)
	X = channel.factorize(A, B, method=method)

	assert_proper_chan(X, m, n)
	assert_arr_equal(A, B @ X, md=1e-4, mrd=0)		# project_subgrad has tolerance 1e-4

def test_factorize_invalid_method(eT, t):
	with pytest.raises(RuntimeError, match="invalid method"):
		channel.factorize(t.id_4, t.noint_4, method="foo")

def test_factorize_subgrad(real_eT, rt):

	# a "fat" B matrix with cols = 1.5 * rows seems to be a good case where the initial solution of A = B X is
	# not a proper matrix, and needs to be improved by the subgradient method (when cols == rows or cols = 2 * rows it
	# seems that the solution is goot right away, with almost no iterations).
	#
	n, m = 10, 15
	B = channel.randu(n, m)
	A = B @ channel.randu(m, n)
	X = channel.factorize(A, B, method="project_subgrad")

	assert_proper_chan(X, m, n)
	assert_arr_equal(A, B @ X, md=1e-4, mrd=0)

	# the following matrices cause the S matrix of the subgradient method to contain inf, causing X to contain -nan
	#
	A = arr(real_eT, "0.4405 0.5595; 0.6588 0.3412")
	B = arr(real_eT, "0.9694 0.0062 0.0244; 0.5312 0.1401 0.3287")
	X = channel.factorize(A, B, method="project_subgrad")

	assert_proper_chan(X, 3, 2)
	assert_arr_equal(A, B @ X, md=1e-4, mrd=0)

def test_factorize_invalid_method(eT, t):
	with pytest.raises(RuntimeError, match="invalid method: foo"):
		channel.factorize(t.id_4, t.id_4, method="foo")
	if eT is rat:
		with pytest.raises(RuntimeError, match="subgrad method not available for rat"):
			channel.factorize(t.id_4, t.id_4, method="subgrad")

@pytest.mark.parametrize("method", ["auto", "solve_lp", "project_subgrad", "project_lp", "project_qp"])
def test_left_factorize(eT, t, method):
	if method in ("project_subgrad", "project_qp") and eT is rat:
		pytest.skip(method + " not available for rat")

	# non factorizable
	assert_proper_chan(channel.left_factorize(t.id_10, t.noint_10, method=method), 0, 0)
	assert_proper_chan(channel.left_factorize(t.id_4,  t.noint_10, method=method), 0, 0)

	n, m = 4, 6
	B = channel.randu(n, m)
	A = channel.randu(m, n) @ B
	X = channel.left_factorize(A, B, method=method)

	assert_proper_chan(X, m, n)
	assert_arr_equal(A, X @ B, md=1e-4, mrd=0)		# subgrad has tolerance 1e-4

def test_left_factorize_large(real_eT):
	# for >= 1000 elements "auto" uses subgrad, which used to fail for col_stoch (used by left_factorize)
	n, m = 30, 40
	B = channel.randu(n, m)
	A = channel.randu(m, n) @ B
	X = channel.left_factorize(A, B)

	assert_proper_chan(X, m, n)
	assert_arr_equal(A, X @ B, md=1e-4, mrd=0)

def test_hyper(eT):
	# cols 1,2 give the same posterior, col 3 has zero probability
	C = arr(eT, "1 0 0 0; 0 1/2 1/2 0; 0 1/2 1/2 0")
	outer, inners = channel.hyper(C, arr(eT, "1/3 1/3 1/3"))

	assert outer.shape == (2,)
	assert inners.shape[1] == 2
	assert_arr_equal(outer, arr(eT, "1/3 2/3"))

	# posteriors without a prior (uniform)
	P = channel.posteriors(arr(eT, "0.8 0.2; 0.3 0.7"))
	assert type(P[0, 0]) is eT		# rat args should not fall back to the double version
	assert_arr_equal(P, arr(eT, "8/11 2/9; 3/11 7/9"))

	R = channel.reduced(C)
	assert_arr_equal(R, arr(eT, "1 0; 0 1; 0 1"))

def test_bayesian_update(eT, t):
	# the identity channel should produce the real prior in 2 iterations
	pi1, iter1 = channel.iterative_bayesian_update(t.id_10, t.prand_10 @ t.id_10)
	assert type(pi1[0]) is eT		# rat args should not fall back to the double version
	assert iter1 == 2
	assert_arr_equal(pi1, t.prand_10)

	# a non interfering channel should produce the uniform prior in 1 iteration
	pi2, iter2 = channel.iterative_bayesian_update(t.noint_10, t.prand_10 @ t.noint_10)
	assert iter2 == 1
	assert_arr_equal(pi2, t.unif_10)

	if eT is double:
		# the geometric should produce the real prior in many iterations (and with limited accuracy)
		# Note: for rat this is slow (probably has to do with the huge denominators in the random elements)
		C = mechanism.d_privacy.geometric(10)
		pi3 = channel.iterative_bayesian_update(C, t.prand_10 @ C, max_diff=1e-8)[0]
		assert_arr_equal(pi3, t.prand_10, md=0, mrd=1e-4)
