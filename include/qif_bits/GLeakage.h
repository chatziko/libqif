/*! \class GLeakage
 *  \brief The Generalized Gain Function model of entropy.
 *
 *  For most information about the foundations of this theory see <a href="../papers/gleakage.pdf">here</a>
 */
template<typename eT>
class GLeakage : public LeakageMeasure<eT> {
	public:
		Mat<eT> G;

		//! A normal constructor taking 2 arguments.
		/*!
		*/
		using LeakageMeasure<eT>::LeakageMeasure;

		GLeakage<eT>() {}
		GLeakage<eT>(const Chan<eT>& C, const Mat<eT>& G) : LeakageMeasure<eT>(C), G(G) {};

		eT vulnerability(const Prob<eT>& pi);
		eT cond_vulnerability(const Prob<eT>& pi);
		eT bayes_risk(const Prob<eT>& pi);
		arma::ucolvec strategy(const Prob<eT>& pi) const;
		arma::ucolvec bayes_strategy(const Prob<eT>& pi) const;

		eT additive_leakage(const Prob<eT>& pi);

		eT entropy(const Prob<eT>& pi)		{ return -internal::real_ops<eT>::log2(vulnerability(pi));		}
		eT cond_entropy(const Prob<eT>& pi)	{ return -internal::real_ops<eT>::log2(cond_vulnerability(pi));	}

//		void * compare_over_prior(chan& other_channel);
//		void * compare_over_gain(chan& other_channel,Prob<eT>& prior);

		virtual const char* class_name() {
			return "GLeakage";
		}

	protected:
		virtual void check_prior(const Prob<eT>& pi, bool ignore_c = false) const {
			if((this->C.n_rows != pi.n_cols && !ignore_c) || G.n_cols != pi.n_cols)
				throw std::runtime_error("invalid prior size");
		}

};


// max_w sum_x pi[x] G[w, x]
//
template<typename eT>
eT GLeakage<eT>::vulnerability(const Prob<eT>& pi) {
	this->check_prior(pi, true);

	return arma::max(this->G * trans(pi));
}

// sum_y max_w sum_x pi[x] C[x, y] G[w, x]
//
template<typename eT>
eT GLeakage<eT>::cond_vulnerability(const Prob<eT>& pi) {
	this->check_prior(pi);

	eT s = eT(0);
	for(uint y = 0; y < this->C.n_cols; y++)
		s += arma::max(this->G * (trans(pi) % this->C.col(y)));
	return s;
}

template<typename eT>
arma::ucolvec GLeakage<eT>::strategy(const Prob<eT>& pi) const {
	this->check_prior(pi);

	arma::ucolvec strategy(pi.n_elem);
	for(uint y = 0; y < this->C.n_cols; y++)
		(this->G * (trans(pi) % this->C.col(y))).max( strategy.at(y) );

	return strategy;
}

// same as cond_vulnerability but with min. G is assumed to be a loss function
//
template<typename eT>
eT GLeakage<eT>::bayes_risk(const Prob<eT>& pi) {
	this->check_prior(pi);

	eT s = eT(0);
	for(uint y = 0; y < this->C.n_cols; y++)
		s += arma::min(this->G * (trans(pi) % this->C.col(y)));
	return s;
}

template<typename eT>
arma::ucolvec GLeakage<eT>::bayes_strategy(const Prob<eT>& pi) const {
	this->check_prior(pi);

	arma::ucolvec strategy(pi.n_elem);
	for(uint y = 0; y < this->C.n_cols; y++)
		(this->G * (trans(pi) % this->C.col(y))).min( strategy.at(y) );

	return strategy;
}
