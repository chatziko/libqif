/*! \class Shannon
 *  \brief The shannon model of entropy.
 *
 *  For most information about the foundations of this theory see <a href="../papers/p1.pdf">here</a>
 */

template<typename eT>
class Shannon : public LeakageMeasure<eT> {
	public:
		// inherit the constructors from parent (C++11 feature)
		using LeakageMeasure<eT>::LeakageMeasure;

		eT entropy(const Prob<eT>& pi);

		eT cond_entropy(const Prob<eT>& pi);

		eT capacity();

		virtual const char* class_name() {
			return "Shannon";
		}
};


// H(X) = - sum_x pi[x] log2(pi[x])
//
template<typename eT>
eT Shannon<eT>::entropy(const Prob<eT>& pi) {
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
eT Shannon<eT>::cond_entropy(const Prob<eT>& pi) {
	this->check_prior(pi);

	eT Hyx = 0;
	for(uint x = 0; x < this->C.n_rows; x++)
		Hyx += pi.at(x) * entropy(this->C.row(x));

	return Hyx + entropy(pi) - entropy(pi * this->C);
}

//Blahut-Arimoto Algorithm
//
template<typename eT>
eT Shannon<eT>::capacity() {
	uint m = this->C.n_rows;
	uint n = this->C.n_cols;

	Prob<eT> F(m), Px(m), Py(m);
	uniform(Px);

	while(1) {
		// Py = output dist
		Py = Px * this->C;

		// update F
		for(uint i = 0; i < m; i++) {
			eT s = 0;
			for(uint j = 0; j < n; j++) {
				eT el = this->C.at(i, j);
				s += el > 0 ? el * log(el / Py.at(j)) : 0;		// NOTE: this is e-base log, not 2-base!
			}
			F.at(i) = exp(s);
		}

		// check stop condition
		eT d = dot(F, Px);
		eT IL = qif::log2(d);
		eT IU = qif::log2(max(F));

		if(equal(IU, IL, this->max_diff, this->max_rel_diff))
			return IL;

		// update Px
		Px %= F / d;		// % is element-wise mult
	}
}

