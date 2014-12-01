#ifndef _QIF_tests_aux_h_
#define _QIF_tests_aux_h_

#include "gtest/gtest.h"

#include "types.h"
#include "Prob.h"
#include "Channel.h"
#include "aux.h"


template<typename eT>
void expect_channel(const Mat<eT>& m, const Channel<eT>& c) {
	EXPECT_EQ(m.n_rows, c.n_rows);
	EXPECT_EQ(m.n_rows, c.n_rows);

	for(uint i = 0; i < c.n_rows; i++)
		for(uint j = 0; j < c.n_cols; j++)
			EXPECT_TRUE(equal(m.at(i, j), c.at(i, j)));

	EXPECT_TRUE(is_proper(c));
}

template<typename eT>
void expect_channel(const std::string& s, const Channel<eT>& c) {
	expect_channel(Mat<eT>(s), c);
}

template<typename eT>
void expect_channel(uint rn, uint cn, const Channel<eT>& c) {
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
