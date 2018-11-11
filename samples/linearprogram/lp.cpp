#include "qif"

#include <iostream>
#include <list>

using namespace qif::lp;
using namespace std;

typedef qif::MatrixEntry<double> ME;

int main() {
	/*
	maximize –2x + 5y, subject to: 
	x <= 200 
	x >= 100
	y <= 170 
	y >= 80
	x + y >= 200 
	Solution 650 given for (x,y) = (100,170)
	*/

	// LinearProrgam solves the linear program
	// {min/max} dot(c,x)
	// subject to A x {>=|==|<=} b
	// x >= 0
	//
	LinearProgram<double> lp;
	lp.maximize = true;

	// cost function
	lp.c.set_size(2);
	lp.c(0) = -2;
	lp.c(1) = 5;

	// coefficients
	lp.b.set_size(5);
	lp.b(0) = 200;
	lp.b(1) = 100;
	lp.b(2) = 170;
	lp.b(3) = 80;
	lp.b(4) = 200;

	// sense
	lp.sense.set_size(5);
	lp.sense(0) = lp.sense(2) = '<';
	lp.sense(1) = lp.sense(3) = lp.sense(4) = '>';

	// matrix A is:
	// 1 0
	// 1 0
	// 0 1
	// 0 1
	// 1 1
	// batch insert (non-zero entries):
	std::list<ME> entries;
	entries.push_back(ME(0, 0, 1));	// A(0,0) = 1
	entries.push_back(ME(1, 0, 1));
	entries.push_back(ME(2, 1, 1));
	entries.push_back(ME(3, 1, 1));
	entries.push_back(ME(4, 0, 1));
	entries.push_back(ME(4, 1, 1));
	qif::fill_spmat(lp.A, 5, 2, entries);

	// same thing but much slower
//	lp.A.set_size(5, 2);
//	lp.A(0,0) = 1;
//	lp.A(1,0) = 1;
//	lp.A(2,1) = 1;
//	lp.A(3,1) = 1;
//	lp.A(4,0) = 1;
//	lp.A(4,1) = 1;

	// methods: simplex_primal, simplex_dual, interior
	lp.method = method_t::simplex_primal;
	lp.glp_msg_level = msg_level_t::on;

	bool solved = lp.solve();

	cout << lp.A;
	cout << lp.status;
	cout
		<< "solved: " << solved
		<< "\nstatus: " << lp.status
		<< "\nmethod: " << lp.method
		<< "\nsolution:\n" << lp.x
		<< "\noptimum: " << lp.optimum()
		<< "\n";
}

