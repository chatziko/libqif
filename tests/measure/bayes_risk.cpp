#include "tests_aux.h"

using namespace measure;

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class BayesRiskTest : public BaseTest<eT> {};
template <typename eT>
class BayesRiskTestReals : public BaseTest<eT> {};

TYPED_TEST_SUITE_P(BayesRiskTest);
TYPED_TEST_SUITE_P(BayesRiskTestReals);		// tests that run only on double/float


TYPED_TEST_P(BayesRiskTest, Vulnerability) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(5)/10, bayes_risk::prior(t.unif_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(9)/10, bayes_risk::prior(t.unif_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),    bayes_risk::prior(t.dirac_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(2)/10, bayes_risk::prior(t.pi1));
}

TYPED_TEST_P(BayesRiskTest, Post_vulnerability) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),    bayes_risk::posterior(t.unif_2, t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),    bayes_risk::posterior(t.dirac_2, t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),    bayes_risk::posterior(t.pi1, t.id_2));

	if(!std::is_same<eT, float>::value) { // precision error on floats
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),    bayes_risk::posterior(t.unif_10, t.id_10));
	}
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),    bayes_risk::posterior(t.dirac_10, t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),    bayes_risk::posterior(t.pi2, t.id_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(9)/10, bayes_risk::posterior(t.unif_10, t.noint_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(0),    bayes_risk::posterior(t.dirac_10, t.noint_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, bayes_risk::prior(t.pi2), bayes_risk::posterior(t.pi2, t.noint_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, bayes_risk::prior(t.pi3), bayes_risk::posterior(t.pi3, t.c1));	// no change in vulnerability
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(9)/40, bayes_risk::posterior(t.pi4, t.c1));

	ASSERT_ANY_THROW(bayes_risk::posterior(t.unif_2, t.id_10));
}

TYPED_TEST_P(BayesRiskTest, Mult_capacity) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, infinity<eT>(), bayes_risk::mult_capacity(t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, infinity<eT>(), bayes_risk::mult_capacity(t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1),          bayes_risk::mult_capacity(t.noint_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(2),          bayes_risk::mult_capacity(t.c1));

	EXPECT_PRED_FORMAT2(equal2<eT>, bayes_risk::mult_leakage(t.unif_2, t.crand_2), bayes_risk::mult_capacity(t.crand_2)); // max_mult_leakage is given for 2-secret uniform prior
}


// run the BayesRiskTest test-case for all types
//
REGISTER_TYPED_TEST_SUITE_P(BayesRiskTest, Vulnerability, Post_vulnerability, Mult_capacity);

INSTANTIATE_TYPED_TEST_SUITE_P(BayesRisk, BayesRiskTest, AllTypes);

