#include "tests_aux.h"

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class GainTest : public BaseTest<eT> {};
template <typename eT>
class GainTestReals : public BaseTest<eT> {};

TYPED_TEST_CASE_P(GainTest);
TYPED_TEST_CASE_P(GainTestReals);		// tests that run only on double/float


// TODO: test more gain functions. Currently these tests are copies from bayes.cpp and test only the id gain function


TYPED_TEST_P(GainTest, Vulnerability) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(5)/10, g::vulnerability(t.id_2, t.unif_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1)/10, g::vulnerability(t.id_10, t.unif_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        g::vulnerability(t.id_4, t.dirac_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(8)/10, g::vulnerability(t.id_2, t.pi1));
}

TYPED_TEST_P(GainTest, Cond_vulnerability) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        g::post_vulnerability(t.id_2, t.unif_2, t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        g::post_vulnerability(t.id_2, t.dirac_2, t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        g::post_vulnerability(t.id_2, t.pi1, t.id_2));
	
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        g::post_vulnerability(t.id_10, t.unif_10, t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        g::post_vulnerability(t.id_10, t.dirac_10, t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        g::post_vulnerability(t.id_10, t.pi2, t.id_10));
	
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1)/10, g::post_vulnerability(t.id_10, t.unif_10, t.noint_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        g::post_vulnerability(t.id_10, t.dirac_10, t.noint_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, g::vulnerability(t.id_10, t.pi2), g::post_vulnerability(t.id_10, t.pi2, t.noint_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, g::vulnerability(t.id_2, t.pi3), g::post_vulnerability(t.id_2, t.pi3, t.c1)); // no change in entropy
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(31)/40, g::post_vulnerability(t.id_2, t.pi4, t.c1));

	ASSERT_ANY_THROW(g::post_vulnerability(t.id_10, t.unif_2, t.id_10));
}

TYPED_TEST_P(GainTestReals, Mulg_leakage) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, 1,              g::mulg_leakage(t.id_2, t.unif_2, t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, qif::log2(10),  g::mulg_leakage(t.id_10, t.unif_10, t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0,              g::mulg_leakage(t.id_10, t.unif_10, t.noint_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, qif::log2(1.5), g::mulg_leakage(t.id_2, t.unif_2, t.c1));
}

TYPED_TEST_P(GainTest, Add_capacity) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1)/2, g::add_capacity(t.unif_2,  t.id_2,      true));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,       g::add_capacity(t.unif_2,  t.id_2,      false));

	EXPECT_PRED_FORMAT2(equal2<eT>, 0,       g::add_capacity(t.dirac_2, t.id_2,      true));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0,       g::add_capacity(t.dirac_2, t.id_2,      false));

	EXPECT_PRED_FORMAT2(equal2<eT>, 0,       g::add_capacity(t.unif_4,  t.noint_4,   true));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0,       g::add_capacity(t.unif_4,  t.noint_4,   false));

	ASSERT_ANY_THROW(g::add_capacity(t.unif_2, t.id_10));
}

TYPED_TEST_P(GainTestReals, Leakage_ge) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	arma::Mat<eT> G;
	Chan<eT> R;
	Chan<eT> T = t.crand_10 * t.crand_10;

	// cases when B can leak more than A
	//
	EXPECT_FALSE(g::leakage_ge(t.crand_10, t.id_10, G));
	EXPECT_FALSE(G.empty());
	EXPECT_TRUE(g::post_vulnerability(G, t.unif_10, t.crand_10) < g::post_vulnerability(G, t.unif_10, t.id_10));

	EXPECT_FALSE(g::leakage_ge(t.noint_10, t.crand_10, G));
	EXPECT_FALSE(G.empty());
	EXPECT_TRUE(g::post_vulnerability(G, t.unif_10, t.noint_10) < g::post_vulnerability(G, t.unif_10, t.crand_10));

	EXPECT_FALSE(g::leakage_ge(T, t.crand_10, G));
	EXPECT_FALSE(G.empty());
	EXPECT_TRUE(g::post_vulnerability(G, t.unif_10, T) < g::post_vulnerability(G, t.unif_10, t.crand_10));

	// cases when B refines A (cannot leak more)
	//
	EXPECT_TRUE(g::leakage_ge(t.id_10, t.crand_10, G, R));
	EXPECT_TRUE(G.empty());
	EXPECT_PRED_FORMAT4(chan_equal4<eT>, t.id_10*R, t.crand_10, 1e-3, 0);

	EXPECT_TRUE(g::leakage_ge(t.crand_10, t.noint_10, G, R));
	EXPECT_TRUE(G.empty());
	EXPECT_PRED_FORMAT4(chan_equal4<eT>, t.crand_10*R, t.noint_10, 1e-3, 0);

	EXPECT_TRUE(g::leakage_ge(t.crand_10, T, G, R));
	EXPECT_TRUE(G.empty());
	EXPECT_PRED_FORMAT4(chan_equal4<eT>, t.crand_10*R, T, 1e-3, 0);
}

// run the GainTest test-case for all types, and the GainTestReals only for double/float
//
REGISTER_TYPED_TEST_CASE_P(GainTest, Vulnerability, Cond_vulnerability, Add_capacity);
REGISTER_TYPED_TEST_CASE_P(GainTestReals, Mulg_leakage, Leakage_ge);

INSTANTIATE_TYPED_TEST_CASE_P(Gain, GainTest, AllTypes);
INSTANTIATE_TYPED_TEST_CASE_P(GainReals, GainTestReals, NativeTypes);

