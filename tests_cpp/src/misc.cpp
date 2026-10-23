#include "tests_aux.h"
using namespace qif;



TEST(MiscTest, Equal) {
	double eps = 1e-30;
	double inf = infinity<double>();

	EXPECT_TRUE(  equal(1.0,     1.0     ));
	EXPECT_TRUE(  equal(1.0,     1.0+eps ));
	EXPECT_TRUE(  equal(1.0+eps, 1.0     ));
	EXPECT_TRUE(  equal(1.0,     1.0-eps ));
	EXPECT_TRUE(  equal(1.0-eps, 1.0     ));
	EXPECT_TRUE(  equal(1.0-eps, 1.0+eps ));
	EXPECT_TRUE(  equal(1.0+eps, 1.0-eps ));

	EXPECT_TRUE(  equal(0.0,     0.0     ));
	EXPECT_TRUE(  equal(0.0,     0.0+eps ));
	EXPECT_TRUE(  equal(0.0+eps, 0.0     ));
	EXPECT_TRUE(  equal(0.0,     0.0-eps ));
	EXPECT_TRUE(  equal(0.0-eps, 0.0     ));
	EXPECT_TRUE(  equal(-0.0,    0.0     ));
	EXPECT_TRUE(  equal(0.0,     -0.0    ));
	EXPECT_TRUE(  equal(-0.0,    0.0+eps ));
	EXPECT_TRUE(  equal(0.0+eps, -0.0    ));
	EXPECT_TRUE(  equal(-0.0,    0.0-eps ));
	EXPECT_TRUE(  equal(0.0-eps, -0.0    ));

	EXPECT_TRUE(  equal(5.0/5.0, 1.0     ));
	EXPECT_TRUE(  equal(1.0,     5.0/5.0 ));
	EXPECT_TRUE(  equal(5.0-5.0, 0.0     ));
	EXPECT_TRUE(  equal(0.0,     5.0-5.0 ));
	EXPECT_TRUE(  equal(0.5*2,   1.0     ));
	EXPECT_TRUE(  equal(1.0,     2*0.5   ));
	EXPECT_TRUE(  equal(inf,     inf     ));
	EXPECT_TRUE(  equal(-inf,    -inf    ));

	EXPECT_TRUE(! equal(1.0,     2.0     ));
	EXPECT_TRUE(! equal(2.0,     0.0     ));
	EXPECT_TRUE(! equal(1.0,     inf     ));
	EXPECT_TRUE(! equal(inf,     1.0     ));
	EXPECT_TRUE(! equal(1.0,     -inf    ));
	EXPECT_TRUE(! equal(-inf,     1.0    ));
	EXPECT_TRUE(! equal(inf,     -inf    ));
	EXPECT_TRUE(! equal(-inf,     inf    ));

	float feps = 1e-10;
	float finf = infinity<float>();

	EXPECT_TRUE(  equal(1.0f,      1.0f      ));
	EXPECT_TRUE(  equal(1.0f,      1.0f+feps ));
	EXPECT_TRUE(  equal(1.0f+feps, 1.0f      ));
	EXPECT_TRUE(  equal(1.0f,      1.0f-feps ));
	EXPECT_TRUE(  equal(1.0f-feps, 1.0f      ));
	EXPECT_TRUE(  equal(1.0f-feps, 1.0f+feps ));
	EXPECT_TRUE(  equal(1.0f+feps, 1.0f-feps ));

	EXPECT_TRUE( equal(0.0f,      0.0f      ));
	EXPECT_TRUE( equal(0.0f,      0.0f+feps ));
	EXPECT_TRUE( equal(0.0f+feps, 0.0f      ));
	EXPECT_TRUE( equal(0.0f,      0.0f-feps ));
	EXPECT_TRUE( equal(0.0f-feps, 0.0f      ));
	EXPECT_TRUE( equal(-0.0f,     0.0f      ));
	EXPECT_TRUE( equal(0.0f,      -0.0f     ));
	EXPECT_TRUE( equal(-0.0f,     0.0f+feps ));
	EXPECT_TRUE( equal(0.0f+feps, -0.0f     ));
	EXPECT_TRUE( equal(-0.0f,     0.0f-feps ));
	EXPECT_TRUE( equal(0.0f-feps, -0.0f     ));

	EXPECT_TRUE(  equal(5.0f/5.0f, 1.0f      ));
	EXPECT_TRUE(  equal(1.0f,      5.0f/5.0f ));
	EXPECT_TRUE(  equal(5.0f-5.0f, 0.0f      ));
	EXPECT_TRUE(  equal(0.0f,      5.0f-5.0f ));
	EXPECT_TRUE(  equal(0.5f*2,    1.0f      ));
	EXPECT_TRUE(  equal(1.0f,      2*0.5f    ));
	EXPECT_TRUE(  equal(finf,      finf      ));
	EXPECT_TRUE(  equal(-finf,     -finf     ));

	EXPECT_TRUE(! equal(1.0f,      2.0f      ));
	EXPECT_TRUE(! equal(2.0f,      0.0f      ));
	EXPECT_TRUE(! equal(1.0f,      finf      ));
	EXPECT_TRUE(! equal(finf,      1.0f      ));
	EXPECT_TRUE(! equal(1.0f,      -finf     ));
	EXPECT_TRUE(! equal(-finf,     1.0f      ));
	EXPECT_TRUE(! equal(finf,      -finf     ));
	EXPECT_TRUE(! equal(-finf,     finf      ));
}



TEST(MiscTest, RatMinMax) {
	// armadillo's min/max start from std::numeric_limits<eT>::lowest()/max(), which are 0 for rat unless overridden
	rmat neg = { { rat(-1), rat(-3) }, { rat(-2, 3), rat(-5) } };
	rmat pos = { { rat(1), rat(3) }, { rat(2, 3), rat(5) } };

	EXPECT_EQ(rat(-2, 3), neg.max());
	EXPECT_EQ(rat(2, 3), pos.min());
	EXPECT_EQ(rat(-2, 3), rmat(arma::max(neg))(0));
	EXPECT_EQ(rat(2, 3), rmat(arma::min(pos))(0));
}
