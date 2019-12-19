#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
namespace m = qif::mechanism;


void init_mechanism_d_privacy_module(py::module m) {

	m.def("distance_matrix",	m::d_privacy::distance_matrix<double>, "n_rows"_a, "n_cols"_a, "d"_a);

	m.def("geometric",			m::d_privacy::geometric<double>, "n_rows"_a, "epsilon"_a = 1, "n_cols"_a = 0, "first_x"_a = 0, "first_y"_a = 0);

	m.def("geometric",			m::d_privacy::exponential<double>, "n_rows"_a, "d"_a, "n_cols"_a = 0);

	m.def("randomized_response",m::d_privacy::randomized_response<double>, "n_rows"_a, "epsilon"_a = 1, "n_cols"_a = 0);

	m.def("tight_constraints",	overload<uint,Metric<double,uint>>(m::d_privacy::tight_constraints<double>), "n_rows"_a, "d"_a);

	m.def("exact_distance",		m::d_privacy::exact_distance<double>, "n_rows"_a, "d"_a);

}