
#include "qif"

namespace qif{

template class GLeakage<double>;
template class GLeakage<float>;
template class GLeakage<rat>;

template class LeakageMeasure<double>;
template class LeakageMeasure<float>;
template class LeakageMeasure<rat>;

template class Mechanism<double>;
template class Mechanism<float>;

template class MinEntropy<double>;
template class MinEntropy<float>;
template class MinEntropy<rat>;


//template class Metric<double, uint>;
//template class Metric<float, uint>;
//template class Euclidean<double, uint>;
//template class Euclidean<float, uint>;

template class Shannon<double>;
template class Shannon<float>;

template class LinearProgram<double>;
template class LinearProgram<float>;
template class LinearProgram<rat>;



double Guessing::entropy(const prob& pi) {
	prob sorted = sort(pi);

	double sum = 0;
	for(uint x = 0; x < sorted.n_cols; ++x) {
		sum += x * sorted.at(x);
	}
	return sum;
}

double Guessing::cond_entropy(const prob& pi) {
	check_prior(pi);

	double result = 0;
	for(uint y = 0; y < C.n_cols; y++) {
		//create the vector vy = pi(1)* C[1,y] ... pi(x)* C[x,y]
		prob vy = pi % C.col(y);
		result += entropy(vy);
	}
	return result;
}

}
