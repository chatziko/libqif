#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
using namespace qif::measure;


void init_shannon_module(py::module m) {

	m.def("prior",      	shannon::prior<double>, "pi"_a);

	m.def("posterior",     	shannon::posterior<double>, "pi"_a, "C"_a);

	m.def("add_leakage",   	shannon::add_leakage<double>, "pi"_a, "C"_a);

	m.def("mult_leakage",  	shannon::mult_leakage<double>, "pi"_a, "C"_a);

	m.def("add_capacity",  	shannon::add_capacity<double>, "C"_a, "md"_a = def_md<double>, "mrd"_a = def_mrd<double>);

}