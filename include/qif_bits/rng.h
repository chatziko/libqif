
namespace rng {

void set_seed(uint seed) {

	arma::arma_rng::set_seed(seed);
	std::srand(seed);
}

void set_seed_random() {
	arma::arma_rng::set_seed_random();

	arma::uvec v = arma::randi<arma::uvec>(1);
	std::cout << "set random " << v(0) << "\n";
	std::srand(v(0));
}

}