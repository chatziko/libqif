#ifndef _QIF_tests_aux_h_
#define _QIF_tests_aux_h_

#include "gtest/gtest.h"

#include "types.h"
#include "aux.h"


template<typename eT>
void expect_channel(const Mat<eT>& m, const Channel<eT>& c) {
	EXPECT_EQ(m.n_rows, c.n_rows);
	EXPECT_EQ(m.n_rows, c.n_rows);

	for(uint i = 0; i < c.n_rows; i++)
		for(uint j = 0; j < c.n_rows; j++)
			EXPECT_TRUE(equal(m.at(i, j), c.at(i, j)));

	EXPECT_TRUE(c.is_proper());
}

template<typename eT>
void expect_channel(const std::string& s, const Channel<eT>& c) {
	expect_channel(Mat<eT>(s), c);
}

template<typename eT>
void expect_channel(uint rn, uint cn, const Channel<eT>& c) {
	EXPECT_EQ(rn, c.n_rows);
	EXPECT_EQ(cn, c.n_rows);

	EXPECT_TRUE(c.is_proper());
}


#endif
