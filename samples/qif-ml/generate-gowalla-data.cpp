#include <qif>
#include <vector>
#include <iostream>
#include "getopt.h"
using namespace qif;
using namespace std;


// configuration
string city = "san_francisco";
string filename = "/home/vagabond/Downloads/datasets/gowalla/Gowalla_totalCheckins.txt";

std::vector<double> levels = { 2, 4, 8 };	// eps = log(level)/100
uint width = 20;
uint split_factor = 10;
uint border = 7;
double cell_size = 150;
uint samples = 1e7;

// computed
uint width_out = (width + 2*border) * split_factor;
double cell_size_out = cell_size / split_factor;

// cell-point mapping. The in grid is centered in the out grid. The corner of the out grid is (0,0), the corder of the in grid is corner_in
point corner_in(border * cell_size, border * cell_size);

auto c2p_in = geo::cell_to_point<double>(width, cell_size, corner_in);
auto c2p_out = geo::cell_to_point<double>(width_out, cell_size_out);
auto p2c_out = geo::point_to_cell<double>(width_out, cell_size_out);

auto euclid = metric::euclidean<double,point>();



void generate_arimoto(prob& pi, double eps, ofstream &file, ofstream &file2) {

	// distance between inputs and outputs
	auto d_inout = metric::compose(euclid, c2p_in, c2p_out);

	// generate mechanism
	prob out = probab::uniform(width_out * width_out);
	chan C = shannon::min_distortion(pi, out, eps/2 * d_inout, 1e-5, 1e-5);

	file2 << C;

	// draw
	arma::Mat<uint> drawn = channel::draw(pi, C, samples);
	LargeAvg<double> mean;

	for(uint i = 0; i < samples; i++) {
		uint x_i = drawn(i,0);
		uint z_i = drawn(i,1);
		point x = c2p_in(x_i);
		point z = c2p_out(z_i);

		mean.add(euclid(x, z));

		file << x_i << "," << z_i << "," << z.x << "," << z.y << "\n";
	}

	std::cout << "arimoto sampled utility: " << mean.value() << "\n";
	std::cout << "arimoto exact utility: " << utility::expected_distance(d_inout, pi, C) << "\n";
	std::cout << "arimoto bayes risk: " << (1-bayes_vuln::posterior(pi, C)) << "\n";
}

void generate_laplace(prob& pi, double eps, ofstream& file) {

	auto drawn = probab::draw(pi, samples);
	LargeAvg<double> mean;

	for(uint x_i : drawn) {
		point x = c2p_in(x_i);
		point z = x + mechanism::planar_laplace_draw(eps);

		mean.add(euclid(x, z));

		file << x_i << ",," <<  z.x << "," << z.y << "\n";
	}

	std::cout << "laplace utility: " << mean.value() << "\n";
}

void generate_geometric(prob& pi, double eps, ofstream& file) {

	uint out_of_grid = width_out * width_out;		// extra observation meaning "we got outside"

	// for efficiency, we batch sample n secrets, and n observations drawn from (0,0). Each z is added to the corresponding origin secret
	auto drawn_x = probab::draw(pi, samples);
	auto drawn_z = mechanism::planar_geometric_draw(cell_size_out, eps, samples);

	LargeAvg<double> mean;

	for(uint i = 0; i < samples; i++) {
		uint x_i = drawn_x(i);
		point x = c2p_in(x_i);
		point z = x + drawn_z[i];		// drawn_z's are drawn from (0,0), so add x

		// map z back to output grid
		uint z_i;
		try {
			z_i = min(p2c_out(z), out_of_grid);
		} catch(std::runtime_error e) {
			z_i = out_of_grid;
		}

		mean.add(euclid(x, z));

		file << x_i << "," << z_i << "," <<  z.x << "," << z.y << "\n";
	}

	std::cout << "geometric utility: " << mean.value() << "\n";
}

// dumps the whole dataset in Euclidean coordinates but without discretization
void dump_dataset(const vector<gowalla::Entry>& db) {
	ofstream file;

	// dump dataset
	file.open("dataset");
	for(uint x : gowalla::to_grid(db, qif::locations[city], width, width, cell_size)) {
		point p = c2p_in(x);
		file << x << "," << p.x << "," << p.y << "\n";
	}
	file.close();

	// dump dataset, no discretization
	file.open("dataset-nodisc");
	for(auto p : gowalla::project_dummy(db, qif::locations[city], width * cell_size, width * cell_size)) {
		// project_dummy puts the bottom-left corner at (0,0), but we want the _center_ of the bottom-left cell at corner_in
		p = p + corner_in + point(-cell_size/2, -cell_size/2);

		file << p.x << "," << p.y << "\n";
	}
	file.close();
}


int main(int argc, char** argv) {
	qif::rng::set_seed_random();

	// parse command line opts
	string gowalla_path = "/home/vagabond/Downloads/datasets/gowalla/Gowalla_totalCheckins.txt";

	struct option longopts[] = {
		{ "gowalla",       required_argument, NULL,         'g' },
		{ 0, 0, 0, 0 } // mark end
	};
	int c;
	while ((c = getopt_long(argc, argv, ":f:", longopts, NULL)) != -1) {
		switch (c) {
		 case 'g':
			gowalla_path = optarg;
			break;
		}
	}


	string city = "san_francisco";

	auto db = gowalla::read_dataset(gowalla_path);
	prob pi = gowalla::to_grid_prior(db, qif::locations[city], width, width, cell_size);
	// prob pi = probab::uniform(width * width);

	dump_dataset(db);

	// file to write
	ofstream file;
	ofstream file2;

	// draw from different mechanisms for all levels
	for(double l : levels) {
		double eps = log(l)/100;

		cout << "eps = ln(" << l << ")/100\n";

		file.open("arimoto-" + to_string(l) + ".csv");
		file2.open("arimoto-" + to_string(l) + ".channel");
		generate_arimoto(pi, eps, file, file2);
		file.close();
		file2.close();

		file.open("laplace-" + to_string(l) + ".csv");
		generate_laplace(pi, eps, file);
		file.close();

		file.open("geometric-" + to_string(l) + ".csv");
		generate_geometric(pi, eps, file);
		file.close();
	}
}
