#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
using namespace qif::refinement;


void init_refinement_module(py::module m) {

	// TODO: maybe change the C++ types for uniformity
	m.def("refined_by",
		[](const chan& A, const chan& B, std::string method) {
			if(method == "factorize") {
				return py::cast(refined_by(A, B));
			} else if(method == "project") {
				mat G; chan R;
				bool res = refined_by(A, B, G, R);
				return py::cast(std::tuple(res, G, R));
			} else {
				throw std::runtime_error("invalid method: " + method);
			}
		},
		"A"_a, "B"_a, "method"_a = "factorize"
	);
	m.def("refined_by",
		[](const rchan& A, const rchan& B, std::string method) {
			if(method == "factorize") {
				return py::cast(refined_by(A, B));
			} else if(method == "project") {
				throw std::runtime_error("project method not available for rat (needs quadratic programming)");
			} else {
				throw std::runtime_error("invalid method: " + method);
			}
		},
		"A"_a, "B"_a, "method"_a = "factorize"
	);

	m.def("max_refined_by",  	max_refined_by<double>, "A"_a, "B"_a);
	m.def("max_refined_by",  	max_refined_by<rat>,    "A"_a, "B"_a);

	m.def("priv_refined_by", 	priv_refined_by<double>, "A"_a, "B"_a);

	m.def("add_metric",      	add_metric<double>, "pi"_a, "A"_a, "B"_a);
	m.def("add_metric",      	add_metric<rat>,    "pi"_a, "A"_a, "B"_a);

	m.def("add_metric_bound",	add_metric_bound<double>, "pi"_a, "A"_a, "B"_a);
	m.def("add_metric_bound",	add_metric_bound<rat>,    "pi"_a, "A"_a, "B"_a);

}
