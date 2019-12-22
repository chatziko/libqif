#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
using namespace qif::channel::compose;


void init_channel_compose_module(py::module m) {

	m.def("parallel", 				parallel<double>, "A"_a, "B"_a);
	m.def("parallel",		 		parallel<rat>,    "A"_a, "B"_a);

	m.def("repeated_independent", 	repeated_independent<double>, "C"_a, "n"_a);
	m.def("repeated_independent", 	repeated_independent<rat>,    "C"_a, "n"_a);

}