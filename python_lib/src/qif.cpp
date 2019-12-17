#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;


void init_channel_module(py::module);
void init_probab_module(py::module);
void init_metric_module(py::module);
void init_measure_module(py::module);


py::handle double_c, uint_c, rat_c, point_c;


PYBIND11_MODULE(qif, m) {
	arma::arma_rng::set_seed_random();

	// Init mp++'s pybind11 integration
	mppp_pybind11::init();

	// Types
    py::class_<point>(m, "point")
		.def(py::init<double,double>())
		.def_readwrite("x", &point::x)
		.def_readwrite("y", &point::y)
		.def_static("from_polar", &point::from_polar)
        .def(py::self + py::self)
        .def("__repr__", &point::to_string);

	auto np = py::module::import("numpy");
	m.attr("double") = np.attr("float64");
	m.attr("uint")   = np.attr("uint64");
	m.attr("rat")    = py::module::import("fractions").attr("Fraction");

	// global class references
	double_c = m.attr("double");
	uint_c   = m.attr("uint");
	rat_c    = m.attr("rat");
	point_c  = m.attr("point");

	// initialize modules
	init_channel_module(m.def_submodule("channel", ""));
	init_probab_module (m.def_submodule("probab", ""));
	init_metric_module (m.def_submodule("metric", ""));
	init_measure_module(m.def_submodule("measure", ""));

	// numpy formatter, so that rats are nicely displayed
	std::function fmt = [](py::object x) { return py::str(x); };
	np.attr("set_printoptions")("formatter"_a = py::dict("object"_a = fmt));
}
