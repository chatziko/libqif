#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;


void init_channel_module(py::module);
void init_probab_module(py::module);
void init_metric_module(py::module);
void init_measure_module(py::module);
void init_mechanism_module(py::module);
void init_refinement_module(py::module);
void init_utility_module(py::module);
void init_lp_module(py::module);


py::handle def_c, double_c, uint_c, rat_c, point_c;


PYBIND11_MODULE(_qif, m) {

	py::module np = py::module::import("numpy");

	// use np.random.randint to get a seed. Use int32_t instead of uint, cause numpy uses int32 for some reason on windows!
	uint seed = np.attr("random").attr("randint")(std::numeric_limits<int32_t>::max()).cast<uint>();
	arma::arma_rng::set_seed(seed);

	// Init mp++'s pybind11 integration
	mppp_pybind11::init();

	// Types
    py::class_<point>(m, "point")
		.def(py::init<double,double>())
		.def_readwrite("x", &point::x)
		.def_readwrite("y", &point::y)
		.def_static("from_polar", &point::from_polar)
		.def_static("from_cell", &point::from_cell)
        .def(py::self + py::self)
        .def("__repr__", &point::to_string);

	// global class references
	double_c = np.attr("float64");
	uint_c   = np.attr("uint32");
	rat_c    = py::module::import("fractions").attr("Fraction");
	point_c  = m.attr("point");
	def_c    = double_c;

	m.def("set_default_type", [](py::object t) { def_c = t; });

	// initialize modules
	init_channel_module   (m.def_submodule("channel",   ""));
	init_probab_module    (m.def_submodule("probab",    ""));
	init_metric_module    (m.def_submodule("metric",    ""));
	init_measure_module   (m.def_submodule("measure",   ""));
	init_mechanism_module (m.def_submodule("mechanism", ""));
	init_refinement_module(m.def_submodule("refinement",""));
	init_utility_module   (m.def_submodule("utility",   ""));
	init_lp_module        (m.def_submodule("lp",        ""));

#ifdef QIF_VERSION
    m.attr("__version__") = QIF_VERSION;
#else
    m.attr("__version__") = "dev";
#endif

}
