
namespace measure::shannon {

// H(X) = - sum_x pi[x] log2(pi[x])
//
template<typename eT>
eT prior(const Prob<eT>& pi) {
	eT sum_x = 0;
	for(uint x = 0; x < pi.n_cols; x++) {
		eT el = pi.at(x);
		sum_x -= el > 0 ? el * qif::log2(el) : 0;
	}

	return sum_x;
}

// computes H(X|Y)
//
// we use the formula
//    H(X|Y) = H(Y|X) + H(X) - H(Y)
// since H(Y|X) is easier to compute:
//    H(Y|X) = sum_x pi[x] H(C[x,-])   (entropy of row x)
//
template<typename eT>
eT posterior(const Prob<eT>& pi, const Chan<eT>& C) {
	channel::check_prior_size(pi, C);

	eT Hyx = 0;
	for(uint x = 0; x < C.n_rows; x++)
		Hyx += pi.at(x) * prior<eT>(C.row(x));

	return Hyx + prior<eT>(pi) - prior<eT>(pi * C);
}

template<typename eT>
eT add_leakage(const Prob<eT>& pi, const Chan<eT>& C) {
	return prior(pi) - posterior(pi, C);
}

template<typename eT>
eT mult_leakage(const Prob<eT>& pi, const Chan<eT>& C) {
	return prior(pi) / posterior(pi, C);
}

//Blahut-Arimoto Algorithm
//
template<typename eT>
eT add_capacity(const Chan<eT>& C, Prob<eT>& Px, eT md = def_md<eT>, eT mrd = def_mrd<eT>) {
	uint m = C.n_rows;
	uint n = C.n_cols;

	Prob<eT> F(m), Py(m);
	Px = probab::uniform<eT>(m);

	while(1) {
		// Py = output dist
		Py = Px * C;

		// update F
		for(uint i = 0; i < m; i++) {
			eT s = 0;
			for(uint j = 0; j < n; j++) {
				eT el = C.at(i, j);
				s += el > 0 ? el * log(el / Py.at(j)) : 0;		// NOTE: this is e-base log, not 2-base!
			}
			F.at(i) = exp(s);
		}

		// check stop condition
		eT d = dot(F, Px);
		eT IL = qif::log2(d);
		eT IU = qif::log2(max(F));

		if(equal(IU, IL, md, mrd))
			return IL;

		// update Px
		Px %= F / d;		// % is element-wise mult
	}
}

// same without getting back the prior
template<typename eT>
eT add_capacity(const Chan<eT>& C, eT md = def_md<eT>, eT mrd = def_mrd<eT>) {
	Prob<eT> Px;
	return add_capacity<eT>(C, Px, md, mrd);
}


// Variant of the Arimoto-Blahut algorithm that returns a channel C that has the best mutual information among those with the same
// distortion as C (see Cover and Thomas, Sectin 13.8). The exact distorion of C can be controlled by a parameter lambda _embeded_ in dist.
// Note that the algorithm really iterates on output distributions, not on channels.

template<typename eT>
Chan<eT> min_distortion(const Prob<eT>& pi, Prob<eT>& out, Metric<eT,uint> dist, eT md = def_md<eT>, eT mrd = def_mrd<eT>) {

	uint n_rows = pi.n_elem;
	uint n_cols = out.n_elem;
	Mat<eT> C = mechanism::distance_matrix(n_rows, n_cols, dist);

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

} // mechanism
