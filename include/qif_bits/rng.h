
namespace rng {

inline
void set_seed(uint seed) {

	arma::arma_rng::set_seed(seed);
	std::srand(seed);
}

inline
void set_seed_random() {
	arma::arma_rng::set_seed_random();

	arma::uvec v = arma::randi<arma::uvec>(1);
	std::srand(v(0));
}

}