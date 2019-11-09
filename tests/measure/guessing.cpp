#include "tests_aux.h"

using namespace measure;

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class GuessingTest : public BaseTest<eT> {};

TYPED_TEST_SUITE_P(GuessingTest);


TYPED_TEST_P(GuessingTest, Vulnerability) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(15)/10, guessing::prior(t.unif_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(55)/10, guessing::prior(t.unif_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,         guessing::prior(t.dirac_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(12)/10, guessing::prior(t.pi1));
}

TYPED_TEST_P(GuessingTest, Post_entropy) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        guessing::posterior(t.unif_2, t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        guessing::posterior(t.dirac_2, t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        guessing::posterior(t.pi1, t.id_2));

	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        guessing::posterior(t.unif_10, t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        guessing::posterior(t.dirac_10, t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,        guessing::posterior(t.pi2, t.id_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, eT(55)/10, guessing::posterior(t.unif_10, t.noint_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1,         guessing::posterior(t.dirac_10, t.noint_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, guessing::prior(t.pi2), guessing::posterior(t.pi2, t.noint_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, guessing::prior(t.pi3), guessing::posterior(t.pi3, t.c1));	// no change in entropy
	EXPECT_PRED_FORMAT2(equal2<eT>, eT(1225)/1000, guessing::posterior(t.pi4, t.c1));

	ASSERT_ANY_THROW(guessing::posterior(t.unif_2, t.id_10));
}


// run the GuessingTest test-case for all types
//
REGISTER_TYPED_TEST_SUITE_P(GuessingTest, Vulnerability, Post_entropy);

INSTANTIATE_TYPED_TEST_SUITE_P(Guessing, GuessingTest, AllTypes);

