#include "tests_aux.h"

using namespace measure;

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class BayesTest : public BaseTest<eT> {};
template <typename eT>
class BayesTestReals : public BaseTest<eT> {};

TYPED_TEST_SUITE_P(BayesTest);
TYPED_TEST_SUITE_P(BayesTestReals);		// tests that run only on double/float


TYPED_TEST_P(BayesTest, Vulnerability) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(5)/10, bayes_vuln::prior(t.unif_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1)/10, bayes_vuln::prior(t.unif_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1),    bayes_vuln::prior(t.dirac_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(8)/10, bayes_vuln::prior(t.pi1));
}

TYPED_TEST_P(BayesTest, Post_vulnerability) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1),    bayes_vuln::posterior(t.unif_2, t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1),    bayes_vuln::posterior(t.dirac_2, t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1),    bayes_vuln::posterior(t.pi1, t.id_2));

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1),    bayes_vuln::posterior(t.unif_10, t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1),    bayes_vuln::posterior(t.dirac_10, t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1),    bayes_vuln::posterior(t.pi2, t.id_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1)/10, bayes_vuln::posterior(t.unif_10, t.noint_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1),    bayes_vuln::posterior(t.dirac_10, t.noint_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, bayes_vuln::prior(t.pi2), bayes_vuln::posterior(t.pi2, t.noint_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, bayes_vuln::prior(t.pi3), bayes_vuln::posterior(t.pi3, t.c1));	// no change in vulnerability
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(31)/40, bayes_vuln::posterior(t.pi4, t.c1));

	ASSERT_ANY_THROW(bayes_vuln::posterior(t.unif_2, t.id_10));
}

TYPED_TEST_P(BayesTest, Mult_capacity) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(2),     bayes_vuln::mult_capacity(t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(10),    bayes_vuln::mult_capacity(t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1),     bayes_vuln::mult_capacity(t.noint_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(15)/10, bayes_vuln::mult_capacity(t.c1));

	EXPECT_PRED_FORMAT2(equal2<eT>, bayes_vuln::mult_leakage(t.unif_10, t.crand_10), bayes_vuln::mult_capacity(t.crand_10)); // max_mult_leakage is given for uniform prior
}


TYPED_TEST_P(BayesTestReals, Min_entropy_leakage) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, 1,              bayes_vuln::min_entropy_leakage(t.unif_2, t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, qif::log2(10),  bayes_vuln::min_entropy_leakage(t.unif_10, t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0,              bayes_vuln::min_entropy_leakage(t.unif_10, t.noint_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, qif::log2(1.5), bayes_vuln::min_entropy_leakage(t.unif_2, t.c1));
}

TYPED_TEST_P(BayesTestReals, Mult_capacity_bound_cap) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	eT v = std::is_same<eT, float>::value
		? 627991744
		: 627991708.193414211273193359375;

	EXPECT_PRED_FORMAT2(equal2<eT>, v, bayes_vuln::mult_capacity_bound_cap(t.id_4, 1e6));
}

// run the BayesTest test-case for all types, and the BayesTestReals only for double/float
//
REGISTER_TYPED_TEST_SUITE_P(BayesTest, Vulnerability, Post_vulnerability, Mult_capacity);
REGISTER_TYPED_TEST_SUITE_P(BayesTestReals, Min_entropy_leakage, Mult_capacity_bound_cap);

INSTANTIATE_TYPED_TEST_SUITE_P(Bayes, BayesTest, AllTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(BayesReals, BayesTestReals, NativeTypes);

