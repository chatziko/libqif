#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;
using namespace qif::measure;


void init_l_risk_module(py::module m) {

	m.doc() = R"pbdoc(
		:math:`\ell`-risk (loss-based Bayes risk).
	)pbdoc";

	m.def("prior",				overload<const  chan&,              const  prob&>(l_risk::prior<double>), "G"_a, "pi"_a);
	m.def("prior",				overload<const Metric<double,uint>&,const  prob&>(l_risk::prior<double>), "G"_a, "pi"_a);
	m.def("prior",				overload<const rchan&,              const rprob&>(l_risk::prior<rat>   ), "G"_a, "pi"_a);
	m.def("prior",				overload<const Metric<rat,uint>&,   const rprob&>(l_risk::prior<rat>   ), "G"_a, "pi"_a);

	m.def("posterior",			overload<const  chan&,              const  prob&,const  chan&>(l_risk::posterior<double>), "G"_a, "pi"_a, "C"_a);
	m.def("posterior",			overload<const Metric<double,uint>&,const  prob&,const  chan&>(l_risk::posterior<double>), "g"_a, "pi"_a, "C"_a);
	m.def("posterior",			overload<const rchan&,              const rprob&,const rchan&>(l_risk::posterior<rat>   ), "G"_a, "pi"_a, "C"_a);
	m.def("posterior",			overload<const Metric<rat,uint>&,   const rprob&,const rchan&>(l_risk::posterior<rat>   ), "g"_a, "pi"_a, "C"_a);

	m.def("add_leakage",		overload<const  chan&,              const  prob&,const  chan&>(l_risk::add_leakage<double>), "G"_a, "pi"_a, "C"_a);
	m.def("add_leakage",		overload<const Metric<double,uint>&,const  prob&,const  chan&>(l_risk::add_leakage<double>), "g"_a, "pi"_a, "C"_a);
	m.def("add_leakage",		overload<const rchan&,              const rprob&,const rchan&>(l_risk::add_leakage<rat>   ), "G"_a, "pi"_a, "C"_a);
	m.def("add_leakage",		overload<const Metric<rat,uint>&,   const rprob&,const rchan&>(l_risk::add_leakage<rat>   ), "g"_a, "pi"_a, "C"_a);

	m.def("mult_leakage",		overload<const  chan&,              const  prob&,const  chan&>(l_risk::mult_leakage<double>), "G"_a, "pi"_a, "C"_a);
	m.def("mult_leakage",		overload<const Metric<double,uint>&,const  prob&,const  chan&>(l_risk::mult_leakage<double>), "g"_a, "pi"_a, "C"_a);
	m.def("mult_leakage",		overload<const rchan&,              const rprob&,const rchan&>(l_risk::mult_leakage<rat>   ), "G"_a, "pi"_a, "C"_a);
	m.def("mult_leakage",		overload<const Metric<rat,uint>&,   const rprob&,const rchan&>(l_risk::mult_leakage<rat>   ), "g"_a, "pi"_a, "C"_a);

	m.def("strategy",			overload<const  chan&,              const  prob&,const  chan&>(l_risk::strategy<double>), "G"_a, "pi"_a, "C"_a);
	m.def("strategy",			overload<const Metric<double,uint>&,const  prob&,const  chan&>(l_risk::strategy<double>), "g"_a, "pi"_a, "C"_a);
	m.def("strategy",			overload<const rchan&,              const rprob&,const rchan&>(l_risk::strategy<rat>   ), "G"_a, "pi"_a, "C"_a);
	m.def("strategy",			overload<const Metric<rat,uint>&,   const rprob&,const rchan&>(l_risk::strategy<rat>   ), "g"_a, "pi"_a, "C"_a);

	m.def("add_capacity",		l_risk::add_capacity<double>, "pi"_a, "C"_a, "one_spanning_g"_a = false);
	m.def("add_capacity",		l_risk::add_capacity<rat>,    "pi"_a, "C"_a, "one_spanning_g"_a = false);

	m.def("loss_to_gain",		l_risk::loss_to_gain<double>, "n_secrets"_a, "n_guesses"_a, "l"_a);
	m.def("loss_to_gain",		l_risk::loss_to_gain<rat>,    "n_secrets"_a, "n_guesses"_a, "l"_a);

	m.def("l_zero_one",			[](double_c_t) { return l_risk::l_zero_one<double>; }, "type"_a = def_type);
	m.def("l_zero_one",			[](rat_c_t   ) { return l_risk::l_zero_one<rat>;    }, "type"_a = def_type);

	m.def("L_zero_one",			[](uint n, double_c_t) { return l_risk::L_zero_one<double>(n); }, "n_rows"_a, "type"_a = def_type);
	m.def("L_zero_one",			[](uint n, rat_c_t   ) { return l_risk::L_zero_one<rat>(n);    }, "n_rows"_a, "type"_a = def_type);


}