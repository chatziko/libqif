#include "tests_aux.h"

using namespace channel;


// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename T>
class ChanTest : public BaseTest<T> {};
template <typename T>
class ChanTestReals : public BaseTest<T> {};

TYPED_TEST_SUITE_P(ChanTest);
TYPED_TEST_SUITE_P(ChanTestReals);


TYPED_TEST_P(ChanTest, Construct) {
	typedef TypeParam eT;

	const char* s = "1 0 0; 0 1 0";
	Chan<eT> C = { {eT(1), eT(0), eT(0)}, {eT(0), eT(1), eT(0)} };

	EXPECT_PRED_FORMAT2(chan_equal2<eT>, Chan<eT>(s),              C); // char*
	EXPECT_PRED_FORMAT2(chan_equal2<eT>, Chan<eT>(std::string(s)), C); // std::string
	EXPECT_PRED_FORMAT2(chan_equal2<eT>, Chan<eT>(C),              C); // copy

	// malformed channel
	//
	const char* s2 = "1 2; 3 0";
	Chan<eT> C2(s2);

	EXPECT_ANY_THROW( assert_proper(Chan<eT>(s2));              ); // char*
	EXPECT_ANY_THROW( assert_proper(Chan<eT>(std::string(s2))); ); // std::string
	EXPECT_ANY_THROW( assert_proper(Chan<eT>(C2));              ); // Mat
}

TYPED_TEST_P(ChanTest, Identity) {
	typedef TypeParam eT;

	Chan<eT> C;
	C = identity<eT>(0);
	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, C, 0, 0);

	C = identity<eT>(3);
	EXPECT_PRED_FORMAT2(chan_equal2<eT>, C, Chan<eT>("1 0 0; 0 1 0; 0 0 1"));
}

TYPED_TEST_P(ChanTest, Randu) {
	typedef TypeParam eT;

	Chan<eT> C(200, 200);
	randu(C);
	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, C, 200, 200);

	C = randu<eT>(5);
	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, C, 5, 5);

	C = randu<eT>(4, 6);
	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, C, 4, 6);
}

TYPED_TEST_P(ChanTest, Factorize) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	// non factorizable
	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, factorize(t.id_10, t.noint_10), 0, 0);
	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, factorize(t.id_4,  t.noint_10), 0, 0);

	// TODO: factorize_lp (default for small sizes) is unstable under float, it fails half of the time, we should investigate
	if(std::is_same<eT, float>::value) return;

	int n = 4, m = 6;
	Chan<eT>
		B = channel::randu<eT>(n, m),
		A = B * channel::randu<eT>(m, n),
		X1 = factorize(A, B),
		Z1 = B * X1;

	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, X1, m, n);
	EXPECT_PRED_FORMAT2(chan_equal2<eT>, A, Z1);

	// factorize_lp
	//
	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, factorize_lp(t.id_10, t.noint_10), 0, 0);
	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, factorize_lp(t.id_4,  t.noint_10), 0, 0);

	Chan<eT>
		X2 = factorize_lp(A, B),
		Z2 = B * X2;

	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, X2, m, n);
	EXPECT_PRED_FORMAT2(chan_equal2<eT>, A, Z2);
}

TYPED_TEST_P(ChanTestReals, FactorizeSubgrad) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	// non factorizable
	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, factorize_subgrad(t.id_10, t.noint_10), 0, 0);
	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, factorize_subgrad(t.id_4,  t.noint_10), 0, 0);

	// a "fat" B matrix with cols = 1.5 * rows seems to be a good case where the initial solution of A = B X is
	// not a proper matrix, and needs to be improved by the subgradient method (when cols == rows or cols = 2 * rows it
	// seems that the solution is goot right away, with almost no iterations).
	//
	int n = 10, m = 15;
	Chan<eT>
		B = randu<eT>(n, m),
		A = B * randu<eT>(m, n),
		X = factorize_subgrad(A, B),
		Z = B * X;

	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, X, m, n);
	EXPECT_PRED_FORMAT4(chan_equal4<eT>, A, Z, 1e-4, 0);

	// the following matrices cause the S matrix of the subgradient method to contain inf, causing X to contain -nan
	//
	A = "0.4405 0.5595;"
		"0.6588 0.3412 ";
	B = "0.9694 0.0062 0.0244;"
		"0.5312 0.1401 0.3287 ";

	X = factorize_subgrad(A, B),
	Z = B * X;

	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, X, 3, 2);
	EXPECT_PRED_FORMAT4(chan_equal4<eT>, A, Z, 1e-4, 0);
}

TYPED_TEST_P(ChanTest, LeftFactorize) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	// non factorizable
	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, left_factorize(t.id_10, t.noint_10), 0, 0);
	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, left_factorize(t.id_4,  t.noint_10), 0, 0);

	// TODO: factorize_lp (default for small sizes) is unstable under float, it fails half of the time, we should investigate
	if(std::is_same<eT, float>::value) return;

	int n = 4, m = 6;
	Chan<eT>
		B = channel::randu<eT>(n, m),
		A = channel::randu<eT>(m, n) * B,
		X1 = left_factorize(A, B),
		Z1 = X1 * B;

	EXPECT_PRED_FORMAT3(chan_is_proper_size3<eT>, X1, m, n);
	EXPECT_PRED_FORMAT4(chan_equal4<eT>, A, Z1, eT(1e-4), eT(0));		// default is subgrad method, with tolerance 1e-4
}

TYPED_TEST_P(ChanTest, BayesianUpdate) {
	typedef TypeParam eT;
	BaseTest<eT>& t = *this;

	// the identity channel should produce the real prior in 2 iterations
	auto [pi1, iter1] = iterative_bayesian_update<eT>(t.id_10, t.prand_10 * t.id_10);
	EXPECT_EQ(2u, iter1);
	EXPECT_PRED_FORMAT2(prob_equal2<eT>, pi1, t.prand_10);

	// a non interfering channel should produce the uniform prior in 1 iteration
	auto [pi2, iter2] = iterative_bayesian_update<eT>(t.noint_10, t.prand_10 * t.noint_10);
	EXPECT_EQ(1u, iter2);
	EXPECT_PRED_FORMAT2(prob_equal2<eT>, pi2, t.unif_10);

	if(std::is_same<eT, double>::value) {
		// the geometric should produce the real prior in many iterations (and with limited accuracy)
		// Note: for rat this is slow (probably has to do with the huge denominators in the random elements)
		auto C = mechanism::d_privacy::geometric<eT>(10);
		auto pi3 = iterative_bayesian_update<eT>(C, t.prand_10 * C, {}, eT(1e-8)).first;
		EXPECT_PRED_FORMAT4(prob_equal4<eT>, pi3, t.prand_10, eT(0), eT(1e-4));
	}
}


// run ChanTest for all types, ChanTestReals only for native types
//
REGISTER_TYPED_TEST_SUITE_P(ChanTest, Construct, Identity, Randu, Factorize, LeftFactorize, BayesianUpdate);
REGISTER_TYPED_TEST_SUITE_P(ChanTestReals, FactorizeSubgrad);

INSTANTIATE_TYPED_TEST_SUITE_P(Chan, ChanTest, AllTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(Chan, ChanTestReals, NativeTypes);

