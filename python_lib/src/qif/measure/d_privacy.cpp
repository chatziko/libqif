#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
using namespace qif::measure;


void init_d_privacy_module(py::module m) {

	m.def("prior",      		d_privacy::prior<double>, "pi"_a, "d"_a);

	m.def("is_private",  	 	d_privacy::is_private<double>, "C"_a, "d"_a);

	m.def("smallest_epsilon",	d_privacy::smallest_epsilon<double>, "C"_a, "d"_a);

}