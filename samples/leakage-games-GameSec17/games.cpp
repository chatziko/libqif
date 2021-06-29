#include <qif>
#include <vector>
#include <iostream>
#include <iomanip>      // std::setprecision
#include <bitset>
#include "getopt.h"
using namespace qif;
using namespace std;



chan qif_games_channel(const vector<vector<chan>>& Cs, const prob& delta, uint a) {
	uint n_def = delta.n_elem;
	uint n_rows = Cs[0][0].n_rows;
	uint n_cols = Cs[0][0].n_cols;
	chan C(n_rows, n_cols);
	C.fill(0);
	for(uint d = 0; d < n_def; d++)
		C += delta(d) * Cs[a][d];
	return C;
}

double qif_games_payoff(const prob& pi, const vector<vector<chan>>& Cs, const prob& delta, uint a) {
	return measure::bayes_vuln::posterior(pi, qif_games_channel(Cs, delta, a));
}

double qif_games_payoff(const prob& pi, const vector<vector<chan>>& Cs, const prob& delta) {
	uint n_adv = Cs.size();
	double max = 0;
	for(uint a = 0; a < n_adv; a++)
		max = std::max(max, qif_games_payoff(pi, Cs, delta, a));
	return max;
}

vector<vector<chan>> rand_channels(uint n_adv, uint n_def, uint n_rows, uint n_cols) {
	vector<vector<chan>> Cs(n_adv);

	for(uint a = 0; a < n_adv; a++)
	for(uint d = 0; d < n_def; d++)
		Cs[a].push_back(channel::randu(n_rows, n_cols));

	return Cs;
}


vector<vector<chan>> read_channels(string path_template) {
	uint n_adv = 0;
	uint n_def = 0;

	while(mat().load(path_template + "-" + std::to_string(n_adv) + "-0", arma::raw_ascii, false))
		n_adv++;
	while(mat().load(path_template + "-0-" + std::to_string(n_def), arma::raw_ascii, false))
		n_def++;
	// n_adv = 1;

	vector<vector<chan>> Cs(n_adv);
	for(uint a = 0; a < n_adv; a++) {
		for(uint d = 0; d < n_def; d++) {
			std::ostringstream path;
			path << path_template << "-" << a << "-" << d;

			Cs[a].push_back(chan());
			if(!Cs[a][d].load(path.str())) {
				std::cerr << "cannot load " << path.str() << "\n";
				throw std::runtime_error("cannot load " + path.str());
			}
		}
	}

	return Cs;
}

vector<vector<chan>> read_channels(string path_template, prob& pi) {
	std::ostringstream path;
	path << path_template << "-prior";

	if(!pi.load(path.str()))
		throw std::runtime_error("cannot load " + path.str());

	return read_channels(path_template);
}

void qif_games_file(string path_template) {

	prob pi;
	auto Cs = read_channels(path_template, pi);
	uint n_adv = Cs.size();
	uint n_def = Cs[0].size();

	// probab::uniform<double>(pi);
	// pi = "0.0137 0.0548 0.2191 0.4382 0.0002 0.0002 0.0548 0.2191";		// from POST paper, gives uniform delta
	// pi = "0.15 0.15 0.1 0.1 0.1 0.1 0.15 0.15";
	// pi = "0.15 0.15 0.15 0.15 0.10 0.10 0.10 0.10";
	probab::assert_proper(pi);
	cout <<pi;

	prob delta = probab::randu<double>(n_def);

	chan C = Cs[0][0];
	C.fill(0);
	for(uint d = 0; d < n_def; d++)
		C += 1.0/n_def * Cs[0][d];
	cout << "-- C with unif delta\n";
	for(uint x = 0; x < C.n_rows; x++)
		cout << std::bitset<8>(x) << C.row(x);	
	cout << "-- column maxima with unif delta" << arma::max(C, 0);
	C.fill(0);
	for(uint d = 0; d < n_def; d++)
		C += delta(d) * Cs[0][d];
	cout << "-- column maxima with rand delta" << arma::max(C, 0);

	// prob maxima(C.n_cols);
	// uint n_bits = C.n_cols - 1;
	// for(uint i = 0; i <= n_bits; i++)
	// 	maxima(i) = i == n_bits
	// 		? 1
	// 		: std::max(
	// 			1/(std::tgamma(i+1) * (n_bits - i)),	// tgamma(n) = (n-1)!
	// 			std::tgamma(n_bits - 1 - i + 1) / std::tgamma(n_bits-1 + 1)
	// 		); 
	// cout << "-- maxima from formula          " << maxima;

	mat payoffs(n_def, n_adv);
	for(uint a = 0; a < n_adv; a++) {
		for(uint d = 0; d < n_def; d++)
			payoffs(d,a) = measure::bayes_vuln::posterior(pi, Cs[a][d]);
		// payoffs(a,n_def) = qif_games_payoff(pi, Cs, delta, a);
	}
	cout << "--- payoffs def x adv:\n" << payoffs 
		<< "random delta: " << qif_games_payoff(pi, Cs, delta) << "\n"
		<< "-------\n";


	auto res = games::minmax_hidden_bayes(pi, Cs);
	cout << std::setprecision(8)
		<< "result: " << res.first << "\n"
		<< "result verify: " << qif_games_payoff(pi, Cs, res.second) << "\n"
		<< "result with uni: " << qif_games_payoff(pi, Cs, probab::uniform<double>(n_def)) << "\n"
		<< "def strategy: " << res.second << "\n";
}

void compare_methods(string path_template) {

	prob pi;
	auto Cs = read_channels(path_template, pi);
	uint n_adv = Cs.size();
	uint n_def = Cs[0].size();

	// probab::uniform<double>(pi);
	// pi = "0.0137 0.0548 0.2191 0.4382 0.0002 0.0002 0.0548 0.2191";		// from POST paper, gives uniform delta
	// pi = "0.15 0.15 0.1 0.1 0.1 0.1 0.15 0.15";
	// pi = "0.15 0.15 0.15 0.15 0.10 0.10 0.10 0.10";
	probab::assert_proper(pi);
	cout <<pi;

	prob delta = probab::randu<double>(n_def);

	chan C = Cs[0][0];
	C.fill(0);
	for(uint d = 0; d < n_def; d++)
		C += 1.0/n_def * Cs[0][d];
	cout << "-- C with unif delta\n";
	for(uint x = 0; x < C.n_rows; x++)
		cout << std::bitset<8>(x) << C.row(x);	
	cout << "-- column maxima with unif delta" << arma::max(C, 0);
	C.fill(0);
	for(uint d = 0; d < n_def; d++)
		C += delta(d) * Cs[0][d];
	cout << "-- column maxima with rand delta" << arma::max(C, 0);

	// prob maxima(C.n_cols);
	// uint n_bits = C.n_cols - 1;
	// for(uint i = 0; i <= n_bits; i++)
	// 	maxima(i) = i == n_bits
	// 		? 1
	// 		: std::max(
	// 			1/(std::tgamma(i+1) * (n_bits - i)),	// tgamma(n) = (n-1)!
	// 			std::tgamma(n_bits - 1 - i + 1) / std::tgamma(n_bits-1 + 1)
	// 		); 
	// cout << "-- maxima from formula          " << maxima;

	mat payoffs(n_def, n_adv);
	for(uint a = 0; a < n_adv; a++) {
		for(uint d = 0; d < n_def; d++)
			payoffs(d,a) = measure::bayes_vuln::posterior(pi, Cs[a][d]);
		// payoffs(a,n_def) = qif_games_payoff(pi, Cs, delta, a);
	}
	cout << "--- payoffs def x adv:\n" << payoffs 
		<< "random delta: " << qif_games_payoff(pi, Cs, delta) << "\n"
		<< "-------\n";


	auto res = games::minmax_hidden_bayes_both2(pi, Cs);
	cout << std::setprecision(8)
		<< "result: " << res.first << "\n"
		<< "result verify: " << qif_games_payoff(pi, Cs, res.second) << "\n"
		<< "result with uni: " << qif_games_payoff(pi, Cs, probab::uniform<double>(n_def)) << "\n"
		<< "def strategy: " << res.second << "\n";
}

void dp_game(string path_template) {
	auto Cs = read_channels(path_template);
	// auto Cs = rand_channels(10, 10, 50, 50);
	uint n_adv = Cs.size();
	if(n_adv == 0) {
		cout << "no channels were read\n";
		return;
	}
	uint n_def = Cs[0].size();

	cout << "adv strategies: " << n_adv << "\n";
	cout << "def strategies: " << n_def << "\n";

	auto [res, delta] = games::dp_hidden(Cs);

	cout << "result: " << res << "\n";
	cout << "def strategy: " << delta << "\n";
}

int main(int argc, char** argv) {
	qif::rng::set_seed_random();

	// parse command line opts
	string path = "";

	struct option longopts[] = {
		{ "path",       required_argument, NULL,         'g' },
		{ 0, 0, 0, 0 } // mark end
	};
	int c;
	while ((c = getopt_long(argc, argv, ":f:", longopts, NULL)) != -1) {
		switch (c) {
		 case 'g':
			path = optarg;
			break;
		 default:
			std::cerr << "usage: leakage-games-GameSec17 --path=<path>\n";
			return -1;
		}
	}
	if(path == "") {
		std::cerr << "usage: leakage-games-GameSec17 --path=<path>\n";
		return 0;
	}

	// for GameSec'17
	// qif_games_file(path);

	// for journal
	// compare_methods(path);

	// for dp example
	dp_game(path);
}
