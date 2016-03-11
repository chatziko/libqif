/*! \class Guessing
 *  \brief The guessing model of entropy.
 *
 *  For most information about this theory see
 */
template<typename eT>
class Guessing : public LeakageMeasure<eT> {
	public:
		using LeakageMeasure<eT>::LeakageMeasure;

		eT entropy(const Prob<eT>& pi);

		eT cond_entropy(const Prob<eT>& pi);

		virtual const char* class_name() {
			return "Guessing";
		}
};

template<typename eT>
eT Guessing<eT>::entropy(const Prob<eT>& pi) {
	Prob<eT> sorted = sort(pi);

	double sum = 0;
	for(uint x = 0; x < sorted.n_cols; ++x) {
		sum += x * sorted.at(x);
	}
	return sum;
}

template<typename eT>
eT Guessing<eT>::cond_entropy(const Prob<eT>& pi) {
	this->check_prior(pi);

	double result = 0;
	for(uint y = 0; y < this->C.n_cols; y++) {
		//create the vector vy = pi(1)* C[1,y] ... pi(x)* C[x,y]
		prob vy = pi % this->C.col(y);
		result += entropy(vy);
	}
	return result;
}

