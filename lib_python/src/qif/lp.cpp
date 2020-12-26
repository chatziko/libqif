#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;


void init_lp_module(py::module m) {

	m.doc() = R"pbdoc(
		Linear solver.
	)pbdoc";

	// Types
    py::class_<lp::Defaults>(m, "defaults")
		.def_readwrite_static("presolve",  &lp::Defaults::presolve)
		.def_readwrite_static("msg_level", &lp::Defaults::msg_level)
		.def_readwrite_static("method",    &lp::Defaults::method)
		.def_readwrite_static("solver",    &lp::Defaults::solver);

}