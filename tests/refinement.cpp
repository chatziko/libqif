#include "tests_aux.h"

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class RefinementTest : public BaseTest<eT> {};
template <typename eT>
class RefinementTestReals : public BaseTest<eT> {};

TYPED_TEST_CASE_P(RefinementTest);
TYPED_TEST_CASE_P(RefinementTestReals);		// tests that run only on double/float



TYPED_TEST_P(RefinementTestReals, Refined_by_project) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	arma::Mat<eT> G;
	Chan<eT> R;
	Chan<eT> T = t.crand_10 * t.crand_10;

	// cases when B can leak more than A
	//
	EXPECT_FALSE(refinement::refined_by(t.crand_10, t.id_10, G));
	EXPECT_FALSE(G.empty());
	EXPECT_TRUE(g_vuln::posterior(G, t.unif_10, t.crand_10) < g_vuln::posterior(G, t.unif_10, t.id_10));

	EXPECT_FALSE(refinement::refined_by(t.noint_10, t.crand_10, G));
	EXPECT_FALSE(G.empty());
	EXPECT_TRUE(g_vuln::posterior(G, t.unif_10, t.noint_10) < g_vuln::posterior(G, t.unif_10, t.crand_10));

	EXPECT_FALSE(refinement::refined_by(T, t.crand_10, G));
	EXPECT_FALSE(G.empty());
	EXPECT_TRUE(g_vuln::posterior(G, t.unif_10, T) < g_vuln::posterior(G, t.unif_10, t.crand_10));

	// cases when B refines A (cannot leak more)
	//
	EXPECT_TRUE(refinement::refined_by(t.id_10, t.crand_10, G, R));
	EXPECT_TRUE(G.empty());
	EXPECT_PRED_FORMAT4(chan_equal4<eT>, t.id_10*R, t.crand_10, 1e-3, 0);

	EXPECT_TRUE(refinement::refined_by(t.crand_10, t.noint_10, G, R));
	EXPECT_TRUE(G.empty());
	EXPECT_PRED_FORMAT4(chan_equal4<eT>, t.crand_10*R, t.noint_10, 1e-3, 0);

	EXPECT_TRUE(refinement::refined_by(t.crand_10, T, G, R));
	EXPECT_TRUE(G.empty());
	EXPECT_PRED_FORMAT4(chan_equal4<eT>, t.crand_10*R, T, 1e-3, 0);
}

TYPED_TEST_P(RefinementTest, Add_metric) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	Chan<eT> T = t.crand_10 * t.crand_10;

	// cases when B refines A (cannot leak more)
	//
	EXPECT_PRED_FORMAT4(equal4<eT>, eT(0), refinement::add_metric(t.unif_4, t.id_4, t.noint_4), 1e-5, 0);

	if(!t.is_rat) { // random channels too slow for rat, numbers become huge
	EXPECT_PRED_FORMAT4(equal4<eT>, eT(0), refinement::add_metric(t.unif_10, t.id_10, t.crand_10), 1e-5, 0);
	EXPECT_PRED_FORMAT4(equal4<eT>, eT(0), refinement::add_metric(t.unif_10, t.crand_10, t.noint_10), 1e-5, 0);
	EXPECT_PRED_FORMAT4(equal4<eT>, eT(0), refinement::add_metric(t.unif_10, t.crand_10, T), 1e-5, 0);
	}

	// this was causing simplex without Bland rule to loop!
	//
	EXPECT_PRED_FORMAT4(equal4<eT>, eT(18643)/2220000, refinement::add_metric<eT>(
		format_rat<eT>("62/100 3/100 35/100"),
		format_rat<eT>("1/10 2/5 1/10 2/5; 1/5 1/5 3/10 3/10; 1/2 1/10 1/10 3/10"),
		format_rat<eT>("1/5 11/50 29/50; 1/5 2/5  2/5; 7/20 2/5 1/4")
    ), 0, 1e-5);
}

// run the RefinementTest test-case for all types, and the RefinementTestReals only for double/float
//
REGISTER_TYPED_TEST_CASE_P(RefinementTest, Add_metric);
REGISTER_TYPED_TEST_CASE_P(RefinementTestReals, Refined_by_project);

INSTANTIATE_TYPED_TEST_CASE_P(Refinement, RefinementTest, AllTypes);
INSTANTIATE_TYPED_TEST_CASE_P(RefinementReals, RefinementTestReals, NativeTypes);
