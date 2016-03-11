#include "tests_aux.h"

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename T>
class ShannonTest : public ::testing::Test {};

TYPED_TEST_CASE_P(ShannonTest);



TYPED_TEST_P(ShannonTest, Entropy) {
	typedef TypeParam eT;

	Prob<eT> pi;

	pi = uniform<Prob<eT>>(2);
	EXPECT_PRED2(equal2<eT>, 1, shannon::entropy(pi));

	pi = uniform<Prob<eT>>(10);
	EXPECT_PRED2(equal2<eT>, qif::log2(10.0), shannon::entropy(pi));

	pi = Prob<eT>("1 0 0 0");
	EXPECT_PRED2(equal2<eT>, 0, shannon::entropy(pi));

	pi = Prob<eT>("0.2 0.8");
	EXPECT_PRED2(equal2<eT>, 0.721928094887362, shannon::entropy(pi));
}

TYPED_TEST_P(ShannonTest, Cond_entropy) {
	typedef TypeParam eT;

	Chan<eT> C;

	C = identity<Chan<eT>>(2);
	EXPECT_PRED2(equal2<eT>, 0, shannon::cond_entropy(uniform<Prob<eT>>(2), C));
	EXPECT_PRED2(equal2<eT>, 0, shannon::cond_entropy(dirac<Prob<eT>>(2), C));
	EXPECT_PRED2(equal2<eT>, 0, shannon::cond_entropy(Prob<eT>("0.2 0.8"), C));

	C = identity<Chan<eT>>(10);
	EXPECT_PRED2(equal2<eT>, 0, shannon::cond_entropy(uniform<Prob<eT>>(10), C));
	EXPECT_PRED2(equal2<eT>, 0, shannon::cond_entropy(dirac<Prob<eT>>(10), C));
	EXPECT_PRED2(equal2<eT>, 0, shannon::cond_entropy(Prob<eT>("0.2 0.8 0 0 0 0 0 0 0 0"), C));

	no_interference(C);
	EXPECT_PRED2(equal2<eT>, qif::log2(10.0), shannon::cond_entropy(uniform<Prob<eT>>(10), C));
	EXPECT_PRED2(equal2<eT>, 0, shannon::cond_entropy(dirac<Prob<eT>>(10), C));

	Prob<eT> pi = randu<Prob<eT>>(10);
	EXPECT_PRED2(equal2<eT>, shannon::entropy(pi), shannon::cond_entropy(pi, C));

	C = Chan<eT>("0.8 0.2; 0.3 0.7");
	pi = "0.25 0.75";
	EXPECT_PRED2(equal2<eT>, 0.669020059980807, shannon::cond_entropy(pi, C));

	C = identity<Chan<eT>>(10);
	ASSERT_ANY_THROW(shannon::cond_entropy<eT>(uniform<Prob<eT>>(2), C););
}

TYPED_TEST_P(ShannonTest, Capacity) {
	typedef TypeParam eT;

	Chan<eT> C;

	C = identity<Chan<eT>>(2);
	EXPECT_PRED2(equal2<eT>, 1, shannon::add_capacity(C));

	C = identity<Chan<eT>>(10);
	EXPECT_PRED2(equal2<eT>, qif::log2(10), shannon::add_capacity(C));

	eT md =	std::is_same<eT, float>::value ? 1e-6 : def_max_diff<eT>();		// the accuracy of 0 capacity is not great for floats
	C = no_interference<Chan<eT>>(10);
	EXPECT_PRED4(equal<eT>, 0, shannon::add_capacity(C), md, 0.0);

	C = Chan<eT>("0.8 0.2; 0.3 0.7");
	EXPECT_PRED2(equal2<eT>, 0.19123813831431799, shannon::add_capacity(C));

	// symmetric
	C = Chan<eT>(
		".3 .2 .5;"
		".5 .3 .2;"
		".2 .5 .3;"
	);
	double cap = qif::log2(C.n_cols) - shannon::entropy<eT>(C.row(0));
	EXPECT_PRED2(equal2<eT>, cap, shannon::add_capacity(C));

	// weakly symmetric
	C = Chan<eT>(
		"0.333333333 0.166666667 0.5;"
		"0.333333333 0.5         0.166666667;"
	);
	cap = qif::log2(C.n_cols) - shannon::entropy<eT>(C.row(0));
	EXPECT_PRED2(equal2<eT>, cap, shannon::add_capacity(C));
}


// run the ChanTest test-case for double, float
//
REGISTER_TYPED_TEST_CASE_P(ShannonTest, Entropy, Cond_entropy, Capacity);

INSTANTIATE_TYPED_TEST_CASE_P(Shannon, ShannonTest, NativeTypes);
