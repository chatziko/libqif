#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
using namespace qif::metric::optimize;


void init_metric_optimize_module(py::module m) {

	m.doc() = R"pbdoc(
		Metric optimization problems.
	)pbdoc";

	m.def("l1_diameter",			l1_diameter<double>, "C"_a, "method"_a = "direct");
	m.def("l1_diameter",			l1_diameter<rat>,    "C"_a, "method"_a = "direct");

	m.def("min_l2_enclosing_ball",	min_l2_enclosing_ball<double>, "C"_a);

}