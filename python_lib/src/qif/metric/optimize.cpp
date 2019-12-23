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

	m.def("l1_diameter",					l1_diameter<double>, "C"_a, "method"_a = "direct");
	m.def("l1_diameter",					l1_diameter<rat>,    "C"_a, "method"_a = "direct");

	m.def("l2_min_enclosing_ball",			l2_min_enclosing_ball<double>, "C"_a);

	m.def("simplex_l1_min_enclosing_ball",	simplex_l1_min_enclosing_ball<double>, "C"_a, "method"_a = "lp", "in_conv_hull"_a = false);
	m.def("simplex_l1_min_enclosing_ball",	simplex_l1_min_enclosing_ball<rat>,    "C"_a, "method"_a = "lp", "in_conv_hull"_a = false);

	m.def("simplex_project", 				simplex_project<double>, "pi"_a);
	m.def("simplex_project", 				simplex_project<rat>,    "pi"_a);

}