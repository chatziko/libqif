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

@pytest.mark.parametrize("method", ["auto", "solve_lp", "project_subgrad", "project_qp"])
def test_refined_by(real_eT, rt, method):
	not_refined, refined = _refinement_cases(rt)

	for A, B in not_refined:
		assert not refinement.refined_by(A, B, method)
	for A, B in refined:
		assert refinement.refined_by(A, B, method)

def test_refined_by_rat(eT, t):
	if eT is not rat:
		pytest.skip("rat only")

	assert refinement.refined_by(t.id_4, t.noint_4)
	assert not refinement.refined_by(t.noint_4, t.id_4)
	assert refinement.refined_by(t.id_4, t.noint_4, "solve_lp")

	with pytest.raises(RuntimeError, match="not available for rat"):
		refinement.refined_by(t.id_4, t.noint_4, "project_qp")
	with pytest.raises(RuntimeError, match="available for rat"):
		refinement.refined_by_witness(t.id_4, t.noint_4)

def test_refined_by_witness(real_eT, rt):
	# "auto" uses project_qp, the only method giving a witness
	not_refined, refined = _refinement_cases(rt)

	for A, B in not_refined:
		res, G, R = refinement.refined_by_witness(A, B)
		assert not res
		assert G.size > 0 and G.min() >= 0 and G.max() <= 1
		assert g_vuln.posterior(G, rt.unif_10, A) < g_vuln.posterior(G, rt.unif_10, B)

	for A, B in refined:
		res, G, R = refinement.refined_by_witness(A, B)
		assert res
		assert G.size == 0
		assert_arr_equal(A @ R, B, md=1e-3, mrd=0)

	with pytest.raises(RuntimeError, match="does not give a witness"):
		refinement.refined_by_witness(rt.id_10, rt.crand_10, "solve_lp")

def test_refined_by_deprecated(real_eT, rt):
	with pytest.warns(DeprecationWarning):
		assert refinement.refined_by(rt.id_10, rt.crand_10, "factorize")

	with pytest.warns(DeprecationWarning):
		res, G, R = refinement.refined_by(rt.crand_10, rt.id_10, "project")
	assert not res
	assert G.size > 0

@pytest.mark.parametrize("method", ["auto", "solve_lp", "project_subgrad", "project_qp"])
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
