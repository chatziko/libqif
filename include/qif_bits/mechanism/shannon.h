
namespace mechanism::shannon {

// Variant of the Arimoto-Blahut algorithm that returns a channel C that has the best mutual information among those with the same
// distortion as C (see Cover and Thomas, Sectin 13.8). The exact distorion of C can be controlled by a parameter lambda _embeded_ in dist.
// Note that the algorithm really iterates on output distributions, not on channels.

template<typename eT>
Chan<eT> min_distortion(const Prob<eT>& pi, Prob<eT>& out, Metric<eT,uint> dist, eT md = def_md<eT>, eT mrd = def_mrd<eT>) {

	uint n_rows = pi.n_elem;
	uint n_cols = out.n_elem;
	Mat<eT> C = mechanism::d_priv::distance_matrix(n_rows, n_cols, dist);

	// uint cnt = 1;
	while(true) {
		Prob<eT> new_out(n_cols);
		new_out.fill(0);

		// row-by-row to avoid constructing a new channel each time
		for(uint x = 0; x < n_rows; x++) {
			Prob<eT> row  =  C.row(x) % out;
			row /= arma::accu(row);
			new_out += pi(x) * row;
		}

		bool done = true;
		for(uint y = 0; y < n_cols; y++)
			done = done && equal(out(y), new_out(y), md, mrd);

		out = new_out;
		// std::cout << " " << cnt++; std::cout.flush();
		// cnt++; if(done) std::cout << cnt << " iterations\n";
		if(done) {
			C.each_row() %= out;
			C.each_col() /= arma::sum(C, 1);
			return C;
		}
	}
}

template<typename eT>
void min_distortion_fast(const Prob<eT>& pi, Prob<eT>& out, Metric<eT,uint> dist, eT md = def_md<eT>, eT mrd = def_mrd<eT>) {

	uint n_rows = pi.n_elem;
	uint n_cols = out.n_elem;
	// Mat<eT> C = mechanism::distance_matrix(n_rows, n_cols, dist);
	uint window = 10;
	Prob<eT> temp(2 * window + 1);

	// uint cnt = 1;
	// while(cnt < 1000) {
		Prob<eT> new_out(n_cols);
		new_out.fill(0);

		// row-by-row to avoid constructing a new channel each time
		for(uint x = 0; x < n_rows; x++) {
			temp.fill(0);
			uint z = 0;
			for(uint y = (x > window ? x-window : 0); y < n_cols && y <= x + window; y++) {
				temp(z++) = out(y) * qif::exp(-dist(x,y));
				// if(cnt == 1)
				// 	std::cout << x << ", " << y << ", " << out(y) << ", " << dist(x,y) << ", " << exp(-dist(x,y))  << "\n";
			}
			// if(cnt == 1) std::cout << temp << "\n";
			temp /= arma::accu(temp);
			z = 0;
			for(uint y = (x > window ? x-window : 0); y < n_cols && y <= x + window; y++) {
				new_out(y) += pi(x) * temp(z++);
			}
		}
		// new_out /= arma::accu(new_out);

		bool done = true;
		for(uint y = 0; y < n_cols; y++)
			done = done && equal(out(y), new_out(y), md, mrd);

		out = new_out;
		// cnt++; if(done) std::cout << cnt << " iterations\n";
		if(done) return;
	// }
}

} // mechanism::shannon
