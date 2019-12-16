#include <qif>
#include "pybind11_aux.h"

namespace py = pybind11;
using namespace py::literals;
using namespace qif;


void init_channel_module(py::module m) {

	m.def("normalize",      overload<const  chan&>(channel::normalize<double>), "C"_a);
	m.def("normalize",      overload<const rchan&>(channel::normalize<rat>   ), "C"_a);

	m.def("identity",		[](uint n, double_c_t) { return channel::identity<double>(n); }, "n_rows"_a, "type"_a = double_c_t());
	m.def("identity",		[](uint n, rat_c_t   ) { return channel::identity<rat>   (n); }, "n_rows"_a, "type"_a);

	m.def("no_interference",[](uint n, uint m, double_c_t) { return channel::no_interference<double>(n, m); }, "n_rows"_a, "n_cols"_a = 1, "type"_a = double_c_t());
	m.def("no_interference",[](uint n, uint m, rat_c_t   ) { return channel::no_interference<rat>   (n, m); }, "n_rows"_a, "n_cols"_a = 1, "type"_a);

	m.def("randu",			[](uint n, uint m, double_c_t) { return channel::randu<double>(n, m); }, "n_rows"_a, "n_cols"_a = 0, "type"_a = double_c_t());
	m.def("randu",			[](uint n, uint m, rat_c_t   ) { return channel::randu<rat>   (n, m); }, "n_rows"_a, "n_cols"_a = 0, "type"_a);

	m.def("deterministic",	[](std::function<uint(uint)> f, uint n, uint m, double_c_t) { return channel::deterministic<double>(f, n, m); }, "map"_a, "n_rows"_a, "n_cols"_a, "type"_a = double_c_t());
	m.def("deterministic",	[](std::function<uint(uint)> f, uint n, uint m, rat_c_t   ) { return channel::deterministic<rat   >(f, n, m); }, "map"_a, "n_rows"_a, "n_cols"_a, "type"_a);

	m.def("is_proper",      channel::is_proper<double>, "C"_a, "mrd"_a = def_mrd<double>);
	m.def("is_proper",      channel::is_proper<rat>,    "C"_a, "mrd"_a = rat(0));

	m.def("assert_proper",  channel::assert_proper<double>, "C"_a);
	m.def("assert_proper",  channel::assert_proper<rat>,    "C"_a);

	m.def("posterior",      channel::posterior<double>, "C"_a, "pi"_a, "col"_a);
	m.def("posterior",      channel::posterior<rat>,    "C"_a, "pi"_a, "col"_a);

	m.def("posteriors",     channel::posteriors<double>, "C"_a, "pi"_a =  prob());
	m.def("posteriors",     channel::posteriors<rat>,    "C"_a, "pi"_a = rprob());

	m.def("hyper",     		channel::hyper<double>, "C"_a, "pi"_a);
	m.def("hyper",     		channel::hyper<rat>,    "C"_a, "pi"_a);

	m.def("reduced",   		channel::reduced<double>, "C"_a);
	m.def("reduced",   		channel::reduced<rat>,    "C"_a);

	m.def("iterative_bayesian_update", channel::iterative_bayesian_update<double>, "C"_a, "out"_a, "start"_a =  prob(), "max_diff"_a = 1e-6, "max_iter"_a = 0);
	m.def("iterative_bayesian_update", channel::iterative_bayesian_update<rat>,    "C"_a, "out"_a, "start"_a = rprob(), "max_diff"_a = 1e-6, "max_iter"_a = 0);

	// TODO: add "method" to factorize
	m.def("factorize",   	channel::factorize<double>, "A"_a, "B"_a, "col_stoch"_a = false);
	m.def("factorize",   	channel::factorize<rat>,    "A"_a, "B"_a, "col_stoch"_a = false);

	m.def("left_factorize", channel::left_factorize<double>, "A"_a, "B"_a, "col_stoch"_a = false);
	m.def("left_factorize", channel::left_factorize<rat>,    "A"_a, "B"_a, "col_stoch"_a = false);

	m.def("sum_column_min", channel::sum_column_min<double>, "C"_a);
	m.def("sum_column_min", channel::sum_column_min<rat>,    "C"_a);

	m.def("draw",     		overload<const  chan&, const  prob&>(channel::draw<double>), "C"_a, "pi"_a);
	m.def("draw",     		overload<const rchan&, const rprob&>(channel::draw<rat>),    "C"_a, "pi"_a);

	m.def("draw",     		overload<const  chan&, const  prob&, uint>(channel::draw<double>), "C"_a, "pi"_a, "n_samples"_a);
	m.def("draw",     		overload<const rchan&, const rprob&, uint>(channel::draw<rat>),    "C"_a, "pi"_a, "n_samples"_a);

}