#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;


void init_metric_optimize_module(py::module);


void init_metric_module(py::module m) {

	m.doc() = R"pbdoc(
		Metrics and other distance functions.

		.. autosummary::
			:toctree: _autosummary
			:template: template.rst

			optimize

		|
	)pbdoc";

	init_metric_optimize_module(m.def_submodule("optimize", ""));


	m.def("euclidean",		[](double_c_t) { return metric::euclidean<double,double>(); }, "type"_a = def_type());
	m.def("euclidean",		[](uint_c_t  ) { return metric::euclidean<double,uint>  (); }, "type"_a = def_type());
	m.def("euclidean",		[](rat_c_t   ) { return metric::euclidean<rat,rat>      (); }, "type"_a = def_type());
	m.def("euclidean",		[](point_c_t ) { return metric::euclidean<double,point> (); }, "type"_a = def_type());

	m.def("euclidean_chain",[](double_c_t) { return metric::euclidean_chain<double>(); }, "type"_a = def_type());
	m.def("euclidean_chain",[](uint_c_t  ) { return metric::euclidean_chain<uint>  (); }, "type"_a = def_type());
	m.def("euclidean_chain",[](rat_c_t   ) { return metric::euclidean_chain<rat>   (); }, "type"_a = def_type());

	m.def("discrete",		[](double_c_t) { return metric::discrete<double,double>(); }, "type"_a = def_type());
	m.def("discrete",		[](uint_c_t  ) { return metric::discrete<double,uint>  (); }, "type"_a = def_type());

	m.def("mult_reals",		[](double_c_t) { return metric::mult_reals<double,double>(); }, "type"_a = def_type());

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

}
