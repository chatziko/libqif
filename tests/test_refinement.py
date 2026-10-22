from qif import rat, refinement
from qif.measure import g_vuln

from helpers import arr, assert_arr_equal, assert_equal, num


def test_refined_by(real_eT, rt):
	T = rt.crand_10 @ rt.crand_10

	# Test both methods: "factorize" (based on channel.factorize) and "project" which also returns
	# a gain function G showing that A is not refined by B, or the factorization R.

	# cases when B can leak more than A
	#
	for A, B in [(rt.crand_10, rt.id_10), (rt.noint_10, rt.crand_10), (T, rt.crand_10)]:
		assert not refinement.refined_by(A, B)

		res, G, R = refinement.refined_by(A, B, "project")
		assert not res
		assert G.size > 0
		assert g_vuln.posterior(G, rt.unif_10, A) < g_vuln.posterior(G, rt.unif_10, B)

	# cases when B refines A (cannot leak more)
	#
	for A, B in [(rt.id_10, rt.crand_10), (rt.crand_10, rt.noint_10), (rt.crand_10, T)]:
		assert refinement.refined_by(A, B)

		res, G, R = refinement.refined_by(A, B, "project")
		assert res
		assert G.size == 0
		assert_arr_equal(A @ R, B, md=1e-3, mrd=0)

def test_add_metric(eT, t):
	# cases when B refines A (cannot leak more)
	#
	assert_equal(num(eT, 0), refinement.add_metric(t.unif_4, t.id_4, t.noint_4)[0], md=1e-5, mrd=0)

	if eT is not rat:		# random channels too slow for rat, numbers become huge
		T = t.crand_10 @ t.crand_10
		assert_equal(0, refinement.add_metric(t.unif_10, t.id_10, t.crand_10)[0], md=1e-5, mrd=0)
		assert_equal(0, refinement.add_metric(t.unif_10, t.crand_10, t.noint_10)[0], md=1e-5, mrd=0)
		assert_equal(0, refinement.add_metric(t.unif_10, t.crand_10, T)[0], md=1e-5, mrd=0)

	# this was causing simplex without Bland rule to loop!
	#
	assert_equal(num(eT, "18643/2220000"), refinement.add_metric(
		arr(eT, "62/100 3/100 35/100"),
		arr(eT, "1/10 2/5 1/10 2/5; 1/5 1/5 3/10 3/10; 1/2 1/10 1/10 3/10"),
		arr(eT, "1/5 11/50 29/50; 1/5 2/5  2/5; 7/20 2/5 1/4"),
	)[0], md=0, mrd=1e-5)
