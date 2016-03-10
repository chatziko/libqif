/*! \class MinEntropy
 *  \brief The Min-Entropy model of entropy.
 *
 *  For most information about the foundations of this theory see <a href="../papers/p1.pdf">here</a>
 */

template<typename eT>
class MinEntropy : public LeakageMeasure<eT> {
	public:
		using LeakageMeasure<eT>::LeakageMeasure;

		eT vulnerability(const Prob<eT>& pi);
		eT cond_vulnerability(const Prob<eT>& pi);
		eT max_mult_leakage();
		arma::ucolvec strategy(const Prob<eT>& pi) const;

		eT entropy(const Prob<eT>& pi)		{ return -internal::real_ops<eT>::log2(vulnerability(pi));		}
		eT cond_entropy(const Prob<eT>& pi)	{ return -internal::real_ops<eT>::log2(cond_vulnerability(pi));	}
		eT capacity()						{ return  internal::real_ops<eT>::log2(max_mult_leakage());		}

		virtual const char* class_name() {
			return "MinEntropy";
		}
};


template<typename eT>
eT MinEntropy<eT>::vulnerability(const Prob<eT>& pi) {
	return arma::max(pi);
}

//sum y max x pi(x) C[x,y]
//
template<typename eT>
eT MinEntropy<eT>::cond_vulnerability(const Prob<eT>& pi) {
	this->check_prior(pi);

	eT s = eT(0);
	for(uint y = 0; y < this->C.n_cols; y++)
		s += arma::max(trans(pi) % this->C.col(y));
	return s;
}

template<typename eT>
arma::ucolvec MinEntropy<eT>::strategy(const Prob<eT>& pi) const {
	this->check_prior(pi);

	arma::ucolvec strategy(pi.n_elem);
	for(uint y = 0; y < this->C.n_cols; y++)
		(trans(pi) % this->C.col(y)).max( strategy.at(y) );

	return strategy;
}

// sum of column maxima
//
template<typename eT>
eT MinEntropy<eT>::max_mult_leakage() {
	return arma::accu(arma::max(this->C, 0));
}

