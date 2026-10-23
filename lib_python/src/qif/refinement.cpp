#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
using namespace qif::refinement;


void init_refinement_module(py::module m) {

	m.doc() = R"pbdoc(
		Refinement relations.
	)pbdoc";

	// methods "factorize" and "project" are deprecated (1.5.0), "project" returns the result of refined_by_witness
	auto refined_by_deprecated = [](std::string& method) -> bool {
		if(method != "factorize" && method != "project")
			return false;

		std::string msg = "refined_by: method=\"" + method + "\" is deprecated, use " +
			(method == "factorize" ? "method=\"auto\"" : "refined_by_witness(A, B, method=\"project_qp\")");
		if(PyErr_WarnEx(PyExc_DeprecationWarning, msg.c_str(), 1) < 0)
			throw py::error_already_set();

		bool project = method == "project";
		method = project ? "project_qp" : "auto";
		return project;
	};

	m.def("refined_by",
		[=](const chan& A, const chan& B, std::string method) -> py::object {
			if(refined_by_deprecated(method)) {
				mat G; chan R;
				bool res = refined_by(A, B, G, R, method);
				return py::cast(std::tuple(res, G, R));
			}
			return py::cast(refined_by(A, B, method));
		},
		"A"_a, "B"_a, "method"_a = "auto"
	);
	m.def("refined_by",
		[=](const rchan& A, const rchan& B, std::string method) -> py::object {
			if(refined_by_deprecated(method)) {
				rmat G; rchan R;
				bool res = refined_by(A, B, G, R, method);
				return py::cast(std::tuple(res, G, R));
			}
			return py::cast(refined_by(A, B, method));
		},
		"A"_a, "B"_a, "method"_a = "auto"
	);

	m.def("refined_by_witness",
		[](const chan& A, const chan& B, const std::string& method) {
			mat G; chan R;
			bool res = refined_by(A, B, G, R, method);
			return std::tuple(res, G, R);
		},
		"A"_a, "B"_a, "method"_a = "auto"
	);
	m.def("refined_by_witness",
		[](const rchan& A, const rchan& B, const std::string& method) {
			rmat G; rchan R;
			bool res = refined_by(A, B, G, R, method);
			return std::tuple(res, G, R);
		},
		"A"_a, "B"_a, "method"_a = "auto"
	);

	m.def("max_refined_by",  	max_refined_by<double>, "A"_a, "B"_a, "method"_a = "auto");
	m.def("max_refined_by",  	max_refined_by<rat>,    "A"_a, "B"_a, "method"_a = "auto");

	m.def("priv_refined_by", 	priv_refined_by<double>, "A"_a, "B"_a);

	m.def("add_metric",      	add_metric<double>, "pi"_a, "A"_a, "B"_a);
	m.def("add_metric",      	add_metric<rat>,    "pi"_a, "A"_a, "B"_a);

	m.def("add_metric_bound",	add_metric_bound<double>, "pi"_a, "A"_a, "B"_a);
	m.def("add_metric_bound",	add_metric_bound<rat>,    "pi"_a, "A"_a, "B"_a);

}
