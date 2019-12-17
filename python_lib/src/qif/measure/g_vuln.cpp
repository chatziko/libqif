#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
using namespace qif::measure;


void init_g_vuln_module(py::module m) {

	m.def("prior",				overload<const  chan&,              const  prob&>(g_vuln::prior<double>), "G"_a, "pi"_a);
	m.def("prior",				overload<const Metric<double,uint>&,const  prob&>(g_vuln::prior<double>), "G"_a, "pi"_a);
	m.def("prior",				overload<const rchan&,              const rprob&>(g_vuln::prior<rat>   ), "G"_a, "pi"_a);
	m.def("prior",				overload<const Metric<rat,uint>&,   const rprob&>(g_vuln::prior<rat>   ), "G"_a, "pi"_a);

	m.def("posterior",			overload<const  chan&,              const  prob&,const  chan&>(g_vuln::posterior<double>), "G"_a, "pi"_a, "C"_a);
	m.def("posterior",			overload<const Metric<double,uint>&,const  prob&,const  chan&>(g_vuln::posterior<double>), "g"_a, "pi"_a, "C"_a);
	m.def("posterior",			overload<const rchan&,              const rprob&,const rchan&>(g_vuln::posterior<rat>   ), "G"_a, "pi"_a, "C"_a);
	m.def("posterior",			overload<const Metric<rat,uint>&,   const rprob&,const rchan&>(g_vuln::posterior<rat>   ), "g"_a, "pi"_a, "C"_a);

	m.def("add_leakage",		overload<const  chan&,              const  prob&,const  chan&>(g_vuln::add_leakage<double>), "G"_a, "pi"_a, "C"_a);
	m.def("add_leakage",		overload<const Metric<double,uint>&,const  prob&,const  chan&>(g_vuln::add_leakage<double>), "g"_a, "pi"_a, "C"_a);
	m.def("add_leakage",		overload<const rchan&,              const rprob&,const rchan&>(g_vuln::add_leakage<rat>   ), "G"_a, "pi"_a, "C"_a);
	m.def("add_leakage",		overload<const Metric<rat,uint>&,   const rprob&,const rchan&>(g_vuln::add_leakage<rat>   ), "g"_a, "pi"_a, "C"_a);

	m.def("mult_leakage",		overload<const  chan&,              const  prob&,const  chan&>(g_vuln::mult_leakage<double>), "G"_a, "pi"_a, "C"_a);
	m.def("mult_leakage",		overload<const Metric<double,uint>&,const  prob&,const  chan&>(g_vuln::mult_leakage<double>), "g"_a, "pi"_a, "C"_a);
	m.def("mult_leakage",		overload<const rchan&,              const rprob&,const rchan&>(g_vuln::mult_leakage<rat>   ), "G"_a, "pi"_a, "C"_a);
	m.def("mult_leakage",		overload<const Metric<rat,uint>&,   const rprob&,const rchan&>(g_vuln::mult_leakage<rat>   ), "g"_a, "pi"_a, "C"_a);

	m.def("strategy",			overload<const  chan&,              const  prob&,const  chan&>(g_vuln::strategy<double>), "G"_a, "pi"_a, "C"_a);
	m.def("strategy",			overload<const Metric<double,uint>&,const  prob&,const  chan&>(g_vuln::strategy<double>), "g"_a, "pi"_a, "C"_a);
	m.def("strategy",			overload<const rchan&,              const rprob&,const rchan&>(g_vuln::strategy<rat>   ), "G"_a, "pi"_a, "C"_a);
	m.def("strategy",			overload<const Metric<rat,uint>&,   const rprob&,const rchan&>(g_vuln::strategy<rat>   ), "g"_a, "pi"_a, "C"_a);

	m.def("add_capacity",		g_vuln::add_capacity<double>, "pi"_a, "C"_a, "one_spanning_g"_a = false);
	m.def("add_capacity",		g_vuln::add_capacity<rat>,    "pi"_a, "C"_a, "one_spanning_g"_a = false);

	m.def("g_id",				[](double_c_t) { return g_vuln::g_id<double>; }, "type"_a = double_c_t());
	m.def("g_id",				[](rat_c_t   ) { return g_vuln::g_id<rat>;    }, "type"_a);

	m.def("G_id",				[](uint n, double_c_t) { return g_vuln::G_id<double>(n); }, "n_rows"_a, "type"_a = double_c_t());
	m.def("G_id",				[](uint n, rat_c_t   ) { return g_vuln::G_id<rat>(n);    }, "n_rows"_a, "type"_a);

	m.def("g_add",				g_vuln::g_add<double>, "G1"_a, "G2"_a);
	m.def("g_add",				g_vuln::g_add<rat>,    "G1"_a, "G2"_a);

	m.def("g_from_posterior",	g_vuln::g_from_posterior<double>, "G"_a, "C"_a);
	m.def("g_from_posterior",	g_vuln::g_from_posterior<rat>,    "G"_a, "C"_a);

}