#include "qif"

#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace qif;
using namespace std;

uint honest;
uint corrupted;

typedef double T;


// Matrix of Crowds protocol with _n honest users, _c corrupted
// and pf probability of forwarding
//
Chan<T> crowds_matrix(uint _n, uint _c, T pf) {

	T n(_n);
	T c(_c);
	T m(n + c);

	Chan<T> C(_n, _n+1);

	T alpha = equal(pf, T(1))
		? T(0)
		: (n - n * pf) / (m - n * pf);
	T beta = (c * (m - pf * (n - 1))) / (m * (m - pf * n));
	T gamma = (c * pf) / (m * (m - pf * n));

	for(uint i = 0; i < n; i++) {
		for(uint j = 0; j < n+1; j++) {
			C.at(i, j) = j == n ? alpha : i == j ? beta : gamma;
		}
	}
	channel::assert_proper(C);

	return C;
}

Mat<T> tiger_g(uint n) {
	Mat<T> G(n+1, n);
	G.submat(0, 0, n-1, n-1).eye();			// first n guesses are the identity
	G.row(n).fill(T(1)/2);
	return G;
}

// biased prior, n users, one has probability p, the others share the remaining 1-p
//
Prob<T> biased_prior(uint n, T p) {
	Prob<T> pi(n);
	pi.fill(equal(p, T(1)) ? T(0) : (1-p)/(n-1));
	pi.at(0) = p;
	probab::assert_proper(pi);
	return pi;
}


void meleakage_by_p() {
	Chan<T>
		C1 = crowds_matrix(honest, corrupted, 0),
		C2 = crowds_matrix(honest, corrupted, 0.5),
		C3 = crowds_matrix(honest, corrupted, 1);

	ofstream file;
	file.open("crowds_data/data-min-p.txt");

	for(T p(0); less_than_or_eq(p, T(1)); p += T(1)/1000) {
		Prob<T> pi = biased_prior(honest, p);

		file << p << "   " 
			<< (bayes::mult_leakage(pi, C1)) << "   "
			<< (bayes::mult_leakage(pi, C2)) << "   "
			<< (bayes::mult_leakage(pi, C3)) << "\n";
	}
}

void meleakage_by_pf() {
	Prob<T> pi1 = biased_prior(honest, T(4)/10);
	Prob<T> pi2 = biased_prior(honest, T(3)/10);
	Prob<T> pi3 = biased_prior(honest, T(2)/10);

	ofstream file;
	file.open("crowds_data/data-min-pf.txt");

	for(T pf(0); less_than_or_eq(pf, T(1)); pf += T(1)/1000) {
		Chan<T> C = crowds_matrix(honest, corrupted, pf);

		file << pf << "   " 
			<< (bayes::mult_leakage(pi1, C)) << "   "
			<< (bayes::mult_leakage(pi2, C)) << "   "
			<< (bayes::mult_leakage(pi3, C)) << "\n";
	}
}

void gleakage_by_p() {
	Chan<T>
		C1 = crowds_matrix(honest, corrupted, 0),
		C2 = crowds_matrix(honest, corrupted, 0.5),
		C3 = crowds_matrix(honest, corrupted, 1);

	// tiger
	Mat<T> G = tiger_g(honest);

	ofstream file;
	file.open("crowds_data/data-tiger-p.txt");

	for(T p(0); less_than_or_eq(p, T(1)); p += T(1)/1000) {
		Prob<T> pi = biased_prior(honest, p);

		file << p << "   " 
			<< (g_vuln::mult_leakage(G, pi, C1)) << "   "
			<< (g_vuln::mult_leakage(G, pi, C2)) << "   "
			<< (g_vuln::mult_leakage(G, pi, C3)) << "\n";
	}
}

void gleakage_by_pf() {
	Prob<T> pi1 = biased_prior(honest, T(4)/10);
	Prob<T> pi2 = biased_prior(honest, T(3)/10);
	Prob<T> pi3 = biased_prior(honest, T(2)/10);

	// tiger
	ofstream file;
	file.open("crowds_data/data-tiger-pf.txt");

	for(T pf(0); less_than_or_eq(pf, T(1)); pf += T(1)/1000) {
		Chan<T> C = crowds_matrix(honest, corrupted, pf);
		Mat<T> G = tiger_g(honest);

		file << pf << "   " 
			<< (g_vuln::mult_leakage(G, pi1, C)) << "   "
			<< (g_vuln::mult_leakage(G, pi2, C)) << "   "
			<< (g_vuln::mult_leakage(G, pi3, C)) << "\n";
	}
}

void repeated_runs() {
	double pf = 0.7;

	Chan<T> C = crowds_matrix(honest, corrupted, pf);
	double pstop = 0.4;

	uint n = 1;
	cout << "channel " << C.n_rows << " x " << C.n_cols << "\n";
	cout << "capacity of single run " <<  bayes::mult_capacity(C) << "\n";
//	cout << "real capacity for " << n << " runs: " << bayes::mult_capacity(channel::comp::repeated_independent(C, n)) << "\n";
	cout << "comp bound for " << n << " runs: " << n * bayes::mult_capacity(C) << "\n";
	cout << "cap_b bound for " << n << " runs: " << bayes::mult_capacity_bound_cap(C, n) << "\n";
	cout << "limit bound for pstop " << pstop << ": " << bayes::mult_capacity(C) / pstop << "\n";

	Chan<T>
		C1 = crowds_matrix(honest, corrupted, 0),
		C2 = crowds_matrix(honest, corrupted, 0.5),
		C3 = crowds_matrix(honest, corrupted, 1);

	ofstream file;
	file.open("crowds_data/repeated.txt");

	T cap1 = bayes::mult_capacity(C1),
	  cap2 = bayes::mult_capacity(C2),
	  cap3 = bayes::mult_capacity(C3);

	for(T pstop(0); less_than_or_eq(pstop, T(1)); pstop += T(1)/1000) {

		file << pstop << "   "
			<< (cap1 / pstop) << "   "
			<< (cap2 / pstop) << "   "
			<< (cap3 / pstop) << "\n";
	}
}

int main() {
	arma::arma_rng::set_seed_random();

	if(system("mkdir -p crowds_data")) {}

	honest = 100;
	corrupted = 5;

	repeated_runs();

	if(system("cd crowds_data && gnuplot ../../../samples/crowds/repeated.plt")) {}	// ignore return value without warning

	// Create 4 plots:
	// - Min-entropy leakage as a function of p for various values of pf
	// - Min-entropy leakage as a function of pf for various values of p
	// - same for g-leakage

	honest = 5;
	corrupted = 1;

	meleakage_by_p();
	meleakage_by_pf();
	gleakage_by_p();
	gleakage_by_pf();

	if(system("cd crowds_data && gnuplot ../../../samples/crowds/crowds.plt")) {}	// ignore return value without warning
}

