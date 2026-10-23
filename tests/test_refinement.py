import pytest
from qif import rat, refinement
from qif.measure import g_vuln

from helpers import arr, assert_arr_equal, assert_equal, num


def _refinement_cases(rt):
	T = rt.crand_10 @ rt.crand_10

	# cases when B can leak more than A, and when B refines A (cannot leak more)
	not_refined = [(rt.crand_10, rt.id_10), (rt.noint_10, rt.crand_10), (T, rt.crand_10)]
	refined     = [(rt.id_10, rt.crand_10), (rt.crand_10, rt.noint_10), (rt.crand_10, T)]
	return not_refined, refined

METHODS     = ["auto", "solve_lp", "project_subgrad", "project_lp", "project_qp"]
RAT_METHODS = ["auto", "solve_lp", "project_lp"]
WITNESS_METHODS = ["auto", "project_lp", "project_qp"]

@pytest.mark.parametrize("method", METHODS)
def test_refined_by(eT, t, method):
	if eT is rat and method not in RAT_METHODS:
		with pytest.raises(RuntimeError, match="not available for rat"):
			refinement.refined_by(t.id_4, t.noint_4, method)
		return

	not_refined, refined = _refinement_cases(t)

	for A, B in not_refined:
		assert not refinement.refined_by(A, B, method)
	for A, B in refined:
		assert refinement.refined_by(A, B, method)

@pytest.mark.parametrize("method", WITNESS_METHODS)
def test_refined_by_witness(eT, t, method):
	if eT is rat and method not in RAT_METHODS:
		with pytest.raises(RuntimeError, match="not available for rat"):
			refinement.refined_by_witness(t.id_4, t.noint_4, method)
		return

	not_refined, refined = _refinement_cases(t)

	for A, B in not_refined:
		res, G, R = refinement.refined_by_witness(A, B, method)
		assert not res
		assert R.size == 0
		assert G.size > 0 and G.min() >= 0 and G.max() <= 1
		assert g_vuln.posterior(G, t.unif_10, A) < g_vuln.posterior(G, t.unif_10, B)

	for A, B in refined:
		res, G, R = refinement.refined_by_witness(A, B, method)
		assert res
		assert G.size == 0
		assert_arr_equal(A @ R, B, md=1e-3, mrd=0)		# exact for rat

	for m in ["solve_lp"] + (["project_subgrad"] if eT is not rat else []):
		with pytest.raises(RuntimeError, match="does not give a witness"):
			refinement.refined_by_witness(t.id_10, t.crand_10, m)

def test_refined_by_deprecated(real_eT, rt):
	with pytest.warns(DeprecationWarning):
		assert refinement.refined_by(rt.id_10, rt.crand_10, "factorize")

	with pytest.warns(DeprecationWarning):
		res, G, R = refinement.refined_by(rt.crand_10, rt.id_10, "project")
	assert not res
	assert G.size > 0

@pytest.mark.parametrize("method", METHODS)
def test_max_refined_by(real_eT, rt, method):
	assert refinement.max_refined_by(rt.id_10, rt.crand_10, method)
	assert not refinement.max_refined_by(rt.noint_10, rt.id_10, method)

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
