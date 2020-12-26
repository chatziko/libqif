
namespace mechanism::shannon {

// Variant of the Arimoto-Blahut algorithm that returns a channel C that has the max posterior entropy (i.e. min mutual information)
// among those with at most the same expected loss (aka distortion) as C (the channel we return). See Cover and Thomas, Sectin 13.8.
// The exact distorion of C can be controlled by scaling "loss" by a parameter lambda (here we assume lambda to be already embbeded in loss)
// Note that the algorithm really iterates on output distributions, not on channels.

template<typename eT>
Chan<eT> max_entropy_given_same_loss(const Prob<eT>& pi, Prob<eT> out, Metric<eT,uint> loss, eT md = def_md<eT>, eT mrd = def_mrd<eT>) {

	uint n_rows = pi.n_elem;
	uint n_cols = out.n_elem;
	Mat<eT> C = mechanism::d_privacy::distance_matrix(n_rows, n_cols, loss);

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

// TODO: clean this up. is it really faster? does it work?
template<typename eT>
void max_entropy_for_same_loss_fast(const Prob<eT>& pi, Prob<eT>& out, Metric<eT,uint> loss, eT md = def_md<eT>, eT mrd = def_mrd<eT>) {

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
				temp(z++) = out(y) * qif::exp(-loss(x,y));
				// if(cnt == 1)
				// 	std::cout << x << ", " << y << ", " << out(y) << ", " << loss(x,y) << ", " << exp(-loss(x,y))  << "\n";
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
