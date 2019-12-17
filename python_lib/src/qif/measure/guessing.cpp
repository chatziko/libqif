#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
using namespace qif::measure;


void init_guessing_module(py::module m) {

	m.def("prior",      	guessing::prior<double>, "pi"_a);
	m.def("prior",      	guessing::prior<rat>,    "pi"_a);

	m.def("posterior",     	guessing::posterior<double>, "pi"_a, "C"_a);
	m.def("posterior",     	guessing::posterior<rat>,    "pi"_a, "C"_a);

	m.def("add_leakage",   	guessing::add_leakage<double>, "pi"_a, "C"_a);
	m.def("add_leakage",   	guessing::add_leakage<rat>,    "pi"_a, "C"_a);

	m.def("mult_leakage",  	guessing::mult_leakage<double>, "pi"_a, "C"_a);
	m.def("mult_leakage",  	guessing::mult_leakage<rat>,    "pi"_a, "C"_a);

}