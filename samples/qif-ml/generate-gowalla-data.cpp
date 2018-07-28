#include <qif>
#include <vector>
#include <iostream>
using namespace qif;
using namespace std;


// configuration
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



void generate_arimoto(prob& pi, double eps, ofstream &file) {

	// distance between inputs and outputs
	auto euclid = metric::euclidean<double,point>();
	auto d_inout = metric::compose(euclid, c2p_in, c2p_out);

	// generate mechanism
	prob out = probab::uniform(width_out * width_out);
	chan C = shannon::min_distortion(pi, out, eps/2 * d_inout, 1e-5, 1e-5);

	// draw
	arma::Mat<uint> drawn = channel::draw(pi, C, samples);

	for(uint i = 0; i < samples; i++) {
		uint x_i = drawn(i,0);
		uint z_i = drawn(i,1);
		point z = c2p_out(z_i);

		file << x_i << "," << z_i << "," << z.x << "," << z.y << "\n";
	}
}

void generate_laplace(prob& pi, double eps, ofstream& file) {

	auto drawn = probab::draw(pi, samples);

	for(uint x_i : drawn) {
		point x = c2p_in(x_i);
		point z = mechanism::planar_laplace_draw_cart(x, eps);

		file << x_i << ",," <<  z.x << "," << z.y << "\n";
	}
}

void generate_geometric(prob& pi, double eps, ofstream& file) {

	uint out_of_grid = width_out * width_out;		// extra observation meaning "we got outside"

	// for efficiency, we batch sample n secrets, and n observations drawn from (0,0). Each z is added to the corresponding origin secret
	auto drawn_x = probab::draw(pi, samples);
	auto drawn_z = mechanism::planar_geometric_draw(point(0,0), cell_size_out, eps, samples);

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

		file << x_i << "," << z_i << "," <<  z.x << "," << z.y << "\n";
	}
}


int main() {
	qif::rng::set_seed_random();

	string city = "san_francisco";
	string filename = "/home/vagabond/Downloads/datasets/gowalla/Gowalla_totalCheckins.txt";

	auto db = gowalla::read_dataset(filename);
	prob pi = gowalla::to_grid_prior(db, qif::locations[city], width, width, cell_size);
	// prob pi = probab::uniform(width * width);

	// file to write
	ofstream file;

	for(double l : levels) {
		double eps = log(l)/100;

		file.open("arimoto-" + to_string(l) + ".csv");
		generate_arimoto(pi, eps, file);
		file.close();

		file.open("laplace-" + to_string(l) + ".csv");
		generate_laplace(pi, eps, file);
		file.close();

		file.open("geometric-" + to_string(l) + ".csv");
		generate_geometric(pi, eps, file);
		file.close();
	}
}
