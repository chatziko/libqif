#include "tests_aux.h"

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class BayesTest : public BaseTest<eT> {};
template <typename eT>
class BayesTestReals : public BaseTest<eT> {};

TYPED_TEST_CASE_P(BayesTest);
TYPED_TEST_CASE_P(BayesTestReals);		// tests that run only on double/float


TYPED_TEST_P(BayesTest, Vulnerability) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal2<eT>, eT(5)/10, bayes::vulnerability(t.unif_2));
	EXPECT_PRED2(equal2<eT>, eT(1)/10, bayes::vulnerability(t.unif_10));
	EXPECT_PRED2(equal2<eT>, 1,        bayes::vulnerability(t.dirac_4));
	EXPECT_PRED2(equal2<eT>, eT(8)/10, bayes::vulnerability(t.pi1));
}

TYPED_TEST_P(BayesTest, Post_vulnerability) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal2<eT>, 1,        bayes::post_vulnerability(t.unif_2, t.id_2));
	EXPECT_PRED2(equal2<eT>, 1,        bayes::post_vulnerability(t.dirac_2, t.id_2));
	EXPECT_PRED2(equal2<eT>, 1,        bayes::post_vulnerability(t.pi1, t.id_2));

	EXPECT_PRED2(equal2<eT>, 1,        bayes::post_vulnerability(t.unif_10, t.id_10));
	EXPECT_PRED2(equal2<eT>, 1,        bayes::post_vulnerability(t.dirac_10, t.id_10));
	EXPECT_PRED2(equal2<eT>, 1,        bayes::post_vulnerability(t.pi2, t.id_10));

	EXPECT_PRED2(equal2<eT>, eT(1)/10, bayes::post_vulnerability(t.unif_10, t.noint_10));
	EXPECT_PRED2(equal2<eT>, 1,        bayes::post_vulnerability(t.dirac_10, t.noint_10));

	EXPECT_PRED2(equal2<eT>, bayes::vulnerability(t.pi2), bayes::post_vulnerability(t.pi2, t.noint_10));

	EXPECT_PRED2(equal2<eT>, bayes::vulnerability(t.pi3), bayes::post_vulnerability(t.pi3, t.c1));	// no change in vulnerability
	EXPECT_PRED2(equal2<eT>, eT(31)/40, bayes::post_vulnerability(t.pi4, t.c1));

	ASSERT_ANY_THROW(bayes::post_vulnerability(t.unif_2, t.id_10));
}

TYPED_TEST_P(BayesTest, Mult_capacity) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal2<eT>, 2,         bayes::mult_capacity(t.id_2));
	EXPECT_PRED2(equal2<eT>, 10,        bayes::mult_capacity(t.id_10));
	EXPECT_PRED2(equal2<eT>, 1,         bayes::mult_capacity(t.noint_10));
	EXPECT_PRED2(equal2<eT>, eT(15)/10, bayes::mult_capacity(t.c1));

	EXPECT_PRED2(equal2<eT>, bayes::mult_leakage(t.unif_10, t.crand_10), bayes::mult_capacity(t.crand_10)); // max_mult_leakage is given for uniform prior
}


TYPED_TEST_P(BayesTestReals, Mulg_leakage) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED2(equal2<eT>, 1,              bayes::mulg_leakage(t.unif_2, t.id_2));
	EXPECT_PRED2(equal2<eT>, qif::log2(10),  bayes::mulg_leakage(t.unif_10, t.id_10));
	EXPECT_PRED2(equal2<eT>, 0,              bayes::mulg_leakage(t.unif_10, t.noint_10));
	EXPECT_PRED2(equal2<eT>, qif::log2(1.5), bayes::mulg_leakage(t.unif_2, t.c1));
}

// run the BayesTest test-case for all types, and the BayesTestReals only for double/float
//
REGISTER_TYPED_TEST_CASE_P(BayesTest, Vulnerability, Post_vulnerability, Mult_capacity);
REGISTER_TYPED_TEST_CASE_P(BayesTestReals, Mulg_leakage);

INSTANTIATE_TYPED_TEST_CASE_P(Bayes, BayesTest, AllTypes);
INSTANTIATE_TYPED_TEST_CASE_P(BayesReals, BayesTestReals, NativeTypes);

