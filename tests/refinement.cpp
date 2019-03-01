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
	EXPECT_TRUE(g::post_vulnerability(G, t.unif_10, t.crand_10) < g::post_vulnerability(G, t.unif_10, t.id_10));

	EXPECT_FALSE(refinement::refined_by(t.noint_10, t.crand_10, G));
	EXPECT_FALSE(G.empty());
	EXPECT_TRUE(g::post_vulnerability(G, t.unif_10, t.noint_10) < g::post_vulnerability(G, t.unif_10, t.crand_10));

	EXPECT_FALSE(refinement::refined_by(T, t.crand_10, G));
	EXPECT_FALSE(G.empty());
	EXPECT_TRUE(g::post_vulnerability(G, t.unif_10, T) < g::post_vulnerability(G, t.unif_10, t.crand_10));

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

// run the RefinementTest test-case for all types, and the RefinementTestReals only for double/float
//
// REGISTER_TYPED_TEST_CASE_P(RefinementTest);
REGISTER_TYPED_TEST_CASE_P(RefinementTestReals, Refined_by_project);

// INSTANTIATE_TYPED_TEST_CASE_P(Refinement, RefinementTest, AllTypes);
INSTANTIATE_TYPED_TEST_CASE_P(RefinementReals, RefinementTestReals, NativeTypes);
