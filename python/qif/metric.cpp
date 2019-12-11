#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;



PYBIND11_MODULE(metric, m) {
	// TODO: move centrally?
	arma::arma_rng::set_seed_random();

	// Init mp++'s pybind11 integration
	mppp_pybind11::init();

	m.def("euclidean",		[](dtype_d) { return metric::euclidean<double,double>(); }, "dtype"_a = float64);
	m.def("euclidean",		[](dtype_u) { return metric::euclidean<double,uint>  (); }, "dtype"_a);
	m.def("euclidean",		[](dtype_r) { return metric::euclidean<rat,rat>      (); }, "dtype"_a);

	m.def("euclidean_chain",[](dtype_d) { return metric::euclidean_chain<double>(); }, "dtype"_a = float64);
	m.def("euclidean_chain",[](dtype_u) { return metric::euclidean_chain<uint>  (); }, "dtype"_a);
	m.def("euclidean_chain",[](dtype_r) { return metric::euclidean_chain<rat>   (); }, "dtype"_a);

	m.def("discrete",		[](dtype_d) { return metric::discrete<double,double>(); }, "dtype"_a = float64);
	m.def("discrete",		[](dtype_u) { return metric::discrete<double,uint>  (); });

	m.def("mult_reals",		[](dtype_d) { return metric::mult_reals<double,double>(); }, "dtype"_a = float64);

	// TODO: Maybe do metric manipulation at the python level to treat all types together
	m.def("scale",  				metric::scale<double,double>);
	m.def("scale",  				metric::scale<rat,rat>);
	m.def("scale",  				metric::scale<double,prob>);

	m.def("min",	  				metric::min<double,double>);
	m.def("max",	  				metric::max<double,double>);
	m.def("mirror",	  				metric::mirror<double,double>);
	m.def("threshold", 				metric::threshold<double,double>);
	m.def("threshold_inf",			metric::threshold_inf<double,double>);

	m.def("from_distance_matrix",	metric::from_distance_matrix<double>);
	m.def("to_distance_matrix",		metric::to_distance_matrix<double>);

	m.def("l1",						metric::l1<double,prob>);
	m.def("l2"		,				metric::l2<double,prob>);

	m.def("total_variation",		metric::total_variation<double,prob>);
	m.def("mult_total_variation",	metric::mult_total_variation<double,prob>);
	m.def("convex_separation_quasi",metric::convex_separation_quasi<double,prob>);
	m.def("convex_separation",		metric::convex_separation<double,prob>);
	m.def("kantorovich",			metric::kantorovich<double,prob>);
	m.def("mult_kantorovich",		metric::mult_kantorovich<double,prob>);


	#ifdef QIF_VERSION
		m.attr("__version__") = QIF_VERSION;
	#else
		m.attr("__version__") = "dev";
	#endif
}
