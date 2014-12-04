#ifndef _QIF_tests_aux_h_
#define _QIF_tests_aux_h_

#include "gtest/gtest.h"

#include "types.h"
#include "Prob.h"
#include "Chan.h"
#include "aux.h"


typedef ::testing::Types<double, float, urat> AllTypes;
typedef ::testing::Types<double, float> NativeTypes;
typedef ::testing::Types<urat> RatTypes;

// define a type-parametrized test case (https://code.google.com/p/googletest/wiki/AdvancedGuide)
template <typename eT>
class BaseTest : public ::testing::Test {
	public:
		Prob<eT>
			unif_2   = uniform<Prob<eT>>(2),
			unif_10  = uniform<Prob<eT>>(10),
			unif_100 = uniform<Prob<eT>>(100),
			dirac_2  = dirac<Prob<eT>>(2),
			dirac_4  = dirac<Prob<eT>>(4),
			dirac_10 = dirac<Prob<eT>>(10),
			pi1      = "0.2 0.8",
			pi2      = "0.2 0.8 0 0 0 0 0 0 0 0",
			pi3      = "0.25 0.75",
			pi4      = "0.75 0.25";

		Chan<eT>
			id_2      = identity<Chan<eT>>(2),
			id_4      = identity<Chan<eT>>(4),
			id_10     = identity<Chan<eT>>(10),
			noint_10  = no_interference<Chan<eT>>(10),
			c1        = "0.8 0.2; 0.3 0.7",
			crand_10  = randu<Chan<eT>>(10),
			crand_100 = randu<Chan<eT>>(100);
};


template<typename eT>
void expect_channel(const Mat<eT>& m, const Chan<eT>& c) {
	EXPECT_EQ(m.n_rows, c.n_rows);
	EXPECT_EQ(m.n_rows, c.n_rows);

	for(uint i = 0; i < c.n_rows; i++)
		for(uint j = 0; j < c.n_cols; j++)
			EXPECT_TRUE(equal(m.at(i, j), c.at(i, j)));

	EXPECT_TRUE(is_proper(c));
}

template<typename eT>
void expect_channel(const std::string& s, const Chan<eT>& c) {
	expect_channel(Mat<eT>(s), c);
}

template<typename eT>
void expect_channel(uint rn, uint cn, const Chan<eT>& c) {
	EXPECT_EQ(rn, c.n_rows);
	EXPECT_EQ(cn, c.n_cols);

	EXPECT_TRUE(is_proper(c));
}


template<typename eT>
void expect_prob(const Prob<eT>& m, const Prob<eT>& p) {
	EXPECT_EQ(m.n_cols, p.n_cols);

	for(uint j = 0; j < p.n_cols; j++)
		EXPECT_TRUE(equal(m.at(j), p.at(j)));

	EXPECT_TRUE(is_proper(p));
}

template<typename eT>
void expect_prob(const std::string& s, const Prob<eT>& p) {
	expect_prob(Prob<eT>(s), p);
}

template<typename eT>
void expect_prob(uint cn, const Prob<eT>& p) {
	EXPECT_EQ(cn, p.n_cols);

	EXPECT_TRUE(is_proper(p));
}


#endif
