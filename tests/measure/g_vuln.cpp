#include "tests_aux.h"

using namespace measure;

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class GainTest : public BaseTest<eT> {};

TYPED_TEST_SUITE_P(GainTest);


// TODO: test more gain functions. Currently these tests are copies from bayes_vuln.cpp and test only the id gain function


TYPED_TEST_P(GainTest, Vulnerability) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(5)/10, g_vuln::prior(t.id_2, t.unif_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1)/10, g_vuln::prior(t.id_10, t.unif_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        g_vuln::prior(t.id_4, t.dirac_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(8)/10, g_vuln::prior(t.id_2, t.pi1));
}

TYPED_TEST_P(GainTest, Post_vulnerability) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        g_vuln::posterior(t.id_2, t.unif_2, t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        g_vuln::posterior(t.id_2, t.dirac_2, t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        g_vuln::posterior(t.id_2, t.pi1, t.id_2));
	
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        g_vuln::posterior(t.id_10, t.unif_10, t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        g_vuln::posterior(t.id_10, t.dirac_10, t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        g_vuln::posterior(t.id_10, t.pi2, t.id_10));
	
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1)/10, g_vuln::posterior(t.id_10, t.unif_10, t.noint_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        g_vuln::posterior(t.id_10, t.dirac_10, t.noint_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, g_vuln::prior(t.id_10, t.pi2), g_vuln::posterior(t.id_10, t.pi2, t.noint_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, g_vuln::prior(t.id_2, t.pi3), g_vuln::posterior(t.id_2, t.pi3, t.c1)); // no change in entropy
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(31)/40, g_vuln::posterior(t.id_2, t.pi4, t.c1));

	ASSERT_ANY_THROW(g_vuln::posterior(t.id_10, t.unif_2, t.id_10));
}

TYPED_TEST_P(GainTest, Add_capacity) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1)/2, g_vuln::add_capacity(t.unif_2,  t.id_2,      true));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,       g_vuln::add_capacity(t.unif_2,  t.id_2,      false));

	EXPECT_PRED_FORMAT2(equal2<eT>, 0,       g_vuln::add_capacity(t.dirac_2, t.id_2,      true));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0,       g_vuln::add_capacity(t.dirac_2, t.id_2,      false));

	EXPECT_PRED_FORMAT2(equal2<eT>, 0,       g_vuln::add_capacity(t.unif_4,  t.noint_4,   true));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0,       g_vuln::add_capacity(t.unif_4,  t.noint_4,   false));

	ASSERT_ANY_THROW(g_vuln::add_capacity(t.unif_2, t.id_10));
}

// run the GainTest test-case for all types
//
REGISTER_TYPED_TEST_SUITE_P(GainTest, Vulnerability, Post_vulnerability, Add_capacity);

INSTANTIATE_TYPED_TEST_SUITE_P(Gain, GainTest, AllTypes);

