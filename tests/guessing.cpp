#include "tests_aux.h"

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class GuessingTest : public BaseTest<eT> {};
template <typename eT>
class GuessingTestReals : public BaseTest<eT> {};

TYPED_TEST_CASE_P(GuessingTest);
TYPED_TEST_CASE_P(GuessingTestReals);		// tests that run only on double/float


TYPED_TEST_P(GuessingTest, Vulnerability) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(15)/10, guessing::entropy(t.unif_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(55)/10, guessing::entropy(t.unif_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,         guessing::entropy(t.dirac_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(12)/10, guessing::entropy(t.pi1));
}

TYPED_TEST_P(GuessingTest, Post_entropy) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        guessing::post_entropy(t.unif_2, t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        guessing::post_entropy(t.dirac_2, t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        guessing::post_entropy(t.pi1, t.id_2));

	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        guessing::post_entropy(t.unif_10, t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        guessing::post_entropy(t.dirac_10, t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        guessing::post_entropy(t.pi2, t.id_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(55)/10, guessing::post_entropy(t.unif_10, t.noint_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,         guessing::post_entropy(t.dirac_10, t.noint_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, guessing::entropy(t.pi2), guessing::post_entropy(t.pi2, t.noint_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, guessing::entropy(t.pi3), guessing::post_entropy(t.pi3, t.c1));	// no change in entropy
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1225)/1000, guessing::post_entropy(t.pi4, t.c1));

	ASSERT_ANY_THROW(guessing::post_entropy(t.unif_2, t.id_10));
}


TYPED_TEST_P(GuessingTestReals, Mulg_leakage) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, qif::log2(1.5), guessing::mulg_leakage(t.unif_2, t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, qif::log2(5.5), guessing::mulg_leakage(t.unif_10, t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0,              guessing::mulg_leakage(t.unif_10, t.noint_10));
}

// run the GuessingTest test-case for all types, and the GuessingTestReals only for double/float
//
REGISTER_TYPED_TEST_CASE_P(GuessingTest, Vulnerability, Post_entropy);
REGISTER_TYPED_TEST_CASE_P(GuessingTestReals, Mulg_leakage);

INSTANTIATE_TYPED_TEST_CASE_P(Guessing, GuessingTest, AllTypes);
INSTANTIATE_TYPED_TEST_CASE_P(GuessingReals, GuessingTestReals, NativeTypes);

