#include "tests_aux.h"

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename T>
class ShannonTest : public BaseTest<T> {};

TYPED_TEST_CASE_P(ShannonTest);



TYPED_TEST_P(ShannonTest, Entropy) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, 1, shannon::prior(t.unif_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, qif::log2(10.0), shannon::prior(t.unif_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0, shannon::prior(t.dirac_4));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0.721928094887362, shannon::prior(t.pi1));
}

TYPED_TEST_P(ShannonTest, Cond_entropy) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	EXPECT_PRED_FORMAT2(equal2<eT>, 0, shannon::posterior(t.unif_2, t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0, shannon::posterior(t.dirac_2, t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0, shannon::posterior(t.pi1, t.id_2));

	EXPECT_PRED_FORMAT2(equal2<eT>, 0, shannon::posterior(t.unif_10, t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0, shannon::posterior(t.dirac_10, t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0, shannon::posterior(t.pi2, t.id_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, qif::log2(10.0), shannon::posterior(t.unif_10, t.noint_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, 0, shannon::posterior(t.dirac_10, t.noint_10));

	Prob<eT> pi = probab::randu<eT>(10);
	EXPECT_PRED_FORMAT2(equal2<eT>, shannon::prior(pi), shannon::posterior(pi, t.noint_10));

	EXPECT_PRED_FORMAT2(equal2<eT>, 0.669020059980807, shannon::posterior(t.pi3, t.c1));

	ASSERT_ANY_THROW(shannon::posterior<eT>(t.unif_2, t.id_10););
}

TYPED_TEST_P(ShannonTest, Capacity) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	Prob<eT> pi;

	EXPECT_PRED_FORMAT2(equal2<eT>, 1, shannon::add_capacity(t.id_2));
	EXPECT_PRED_FORMAT2(equal2<eT>, 1, shannon::add_capacity(t.id_2, pi));
	EXPECT_PRED_FORMAT2(prob_equal2<eT>, pi, t.unif_2);

	EXPECT_PRED_FORMAT2(equal2<eT>, qif::log2(10), shannon::add_capacity(t.id_10));
	EXPECT_PRED_FORMAT2(equal2<eT>, qif::log2(10), shannon::add_capacity(t.id_10, pi));
	EXPECT_PRED_FORMAT2(prob_equal2<eT>, pi, t.unif_10);

	eT md =	std::is_same<eT, float>::value ? 1e-6 : def_md<eT>;		// the accuracy of 0 capacity is not great for floats
	EXPECT_PRED_FORMAT4(equal4<eT>, 0, shannon::add_capacity(t.noint_10), md, 0.0);

	EXPECT_PRED_FORMAT2(equal2<eT>, 0.19123813831431799, shannon::add_capacity(t.c1));

	// symmetric
	Chan<eT> C(
		".3 .2 .5;"
		".5 .3 .2;"
		".2 .5 .3;"
	);
	double cap = qif::log2(C.n_cols) - shannon::prior<eT>(C.row(0));
	EXPECT_PRED_FORMAT2(equal2<eT>, cap, shannon::add_capacity(C));

	// weakly symmetric
	C = Chan<eT>(
		"0.333333333 0.166666667 0.5;"
		"0.333333333 0.5         0.166666667;"
	);
	cap = qif::log2(C.n_cols) - shannon::prior<eT>(C.row(0));
	EXPECT_PRED_FORMAT2(equal2<eT>, cap, shannon::add_capacity(C));
}


// run the ChanTest test-case for double, float
//
REGISTER_TYPED_TEST_CASE_P(ShannonTest, Entropy, Cond_entropy, Capacity);

INSTANTIATE_TYPED_TEST_CASE_P(Shannon, ShannonTest, NativeTypes);

