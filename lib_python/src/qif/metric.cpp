#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;


void init_metric_optimize_module(py::module);


// Metrics on secrets (uint): the distances are of the default type (double or rat)
template<typename D, typename R>
py::object metric_on_uint(D d_double, R d_rat) {
	return def_c.is(rat_c) ? py::cast(d_rat) : py::cast(d_double);
}

void init_metric_module(py::module m) {

	init_metric_optimize_module(m.def_submodule("optimize", ""));


	m.def("euclidean",		[](double_c_t) { return metric::euclidean<double,double>(); }, "type"_a = def_type());
	m.def("euclidean",		[](uint_c_t  ) { return metric_on_uint(metric::euclidean<double,uint>(), metric::euclidean<rat,uint>()); }, "type"_a = def_type());
	m.def("euclidean",		[](rat_c_t   ) { return metric::euclidean<rat,rat>      (); }, "type"_a = def_type());
	m.def("euclidean",		[](point_c_t ) { return metric::euclidean<double,point> (); }, "type"_a = def_type());

	m.def("euclidean_chain",[](double_c_t) { return metric::euclidean_chain<double>(); }, "type"_a = def_type());
	m.def("euclidean_chain",[](uint_c_t  ) { return metric::euclidean_chain<uint>  (); }, "type"_a = def_type());
	m.def("euclidean_chain",[](rat_c_t   ) { return metric::euclidean_chain<rat>   (); }, "type"_a = def_type());

	m.def("discrete",		[](double_c_t) { return metric::discrete<double,double>(); }, "type"_a = def_type());
	m.def("discrete",		[](uint_c_t  ) { return metric_on_uint(metric::discrete<double,uint>(), metric::discrete<rat,uint>()); }, "type"_a = def_type());

	m.def("mult_reals",		[](double_c_t) { return metric::mult_reals<double,double>(); }, "type"_a = def_type());

	// scale, min, max, mirror, threshold, threshold_inf are implemented in python (metric/__init__.py)

	m.def("from_distance_matrix",	metric::from_distance_matrix<double>, "M"_a);
	m.def("from_distance_matrix",	metric::from_distance_matrix<rat>,    "M"_a);
	m.def("to_distance_matrix",		[](Metric<double,uint> d, uint n_rows, uint n_cols, double_c_t) { return metric::to_distance_matrix<double>(d, n_rows, n_cols); }, "d"_a, "n_rows"_a, "n_cols"_a = 0, "type"_a = def_type());
	m.def("to_distance_matrix",		[](Metric<rat,uint>    d, uint n_rows, uint n_cols, rat_c_t   ) { return metric::to_distance_matrix<rat>   (d, n_rows, n_cols); }, "d"_a, "n_rows"_a, "n_cols"_a = 0, "type"_a = def_type());

	m.def("manhattan",				metric::manhattan<double,point>);

	m.def("l1",						[](double_c_t) { return metric::l1<double,prob>(); }, "type"_a = def_type());
	m.def("l1",						[](rat_c_t   ) { return metric::l1<rat,rprob>  (); }, "type"_a = def_type());
	m.def("l2",						metric::l2<double,prob>);

	m.def("total_variation",		[](double_c_t) { return metric::total_variation<double,prob>(); }, "type"_a = def_type());
	m.def("total_variation",		[](rat_c_t   ) { return metric::total_variation<rat,rprob>  (); }, "type"_a = def_type());
	m.def("mult_total_variation",	metric::mult_total_variation<double,prob>);
	m.def("convex_separation_quasi",[](double_c_t) { return metric::convex_separation_quasi<double,prob>(); }, "type"_a = def_type());
	m.def("convex_separation_quasi",[](rat_c_t   ) { return metric::convex_separation_quasi<rat,rprob>  (); }, "type"_a = def_type());
	m.def("convex_separation",		[](double_c_t) { return metric::convex_separation<double,prob>(); }, "type"_a = def_type());
	m.def("convex_separation",		[](rat_c_t   ) { return metric::convex_separation<rat,rprob>  (); }, "type"_a = def_type());
	m.def("kantorovich",			[](Metric<double,uint> d, std::string method, double_c_t) { return metric::kantorovich<double,prob>(d, method); }, "d"_a, "method"_a = "auto", "type"_a = def_type());
	m.def("kantorovich",			[](Metric<rat,uint>    d, std::string method, rat_c_t   ) { return metric::kantorovich<rat,rprob>  (d, method); }, "d"_a, "method"_a = "auto", "type"_a = def_type());
	m.def("mult_kantorovich",		metric::mult_kantorovich<double,prob>, "d"_a);

}
