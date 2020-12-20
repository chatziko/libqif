// Adapted from eigen.h

#pragma once

#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <armadillo>
#include <mp++/mp++.hpp>
#include <mp++/extra/pybind11.hpp>
#include <qif>

using qif::uint;

template<typename... Args>
constexpr auto overload = pybind11::overload_cast<Args...>;	// for selecting member of overloaded function

extern pybind11::handle def_c, double_c, uint_c, rat_c, point_c;

struct double_c_t {};		// a type that only accepts double_c
struct uint_c_t {};
struct rat_c_t {};
struct point_c_t {};

// Use a function with a static variable, instead of a global const, to workaround a MSVC issue
// https://github.com/pybind/pybind11/issues/2493
inline pybind11::object def_type() {
	static pybind11::object def_type = pybind11::cast(1);
    return def_type;
}


namespace pybind11::detail {

template <> struct type_caster<double_c_t> {
public:
	PYBIND11_TYPE_CASTER(double_c_t, _("class[double]"));
	bool load(handle src, bool) {
		return (src.is(def_type()) ? def_c : src).is(double_c);
	}
};
template <> struct type_caster<uint_c_t> {
public:
	PYBIND11_TYPE_CASTER(uint_c_t, _("class[uint]"));
	bool load(handle src, bool) {
		return (src.is(def_type()) ? def_c : src).is(uint_c);
	}
};
template <> struct type_caster<rat_c_t> {
public:
	PYBIND11_TYPE_CASTER(rat_c_t, _("class[fraction]"));
	bool load(handle src, bool) {
		return (src.is(def_type()) ? def_c : src).is(rat_c);
	}
};
template <> struct type_caster<point_c_t> {
public:
	PYBIND11_TYPE_CASTER(point_c_t, _("class[point]"));
	bool load(handle src, bool) {
		return (src.is(def_type()) ? def_c : src).is(point_c);
	}
};

///////////////////////////// Mat<native> | Row<native> | Col<native> ///////////////////////////////////////

template <typename Type>
class type_caster<
	Type,
	typename std::enable_if_t<
		arma::is_Mat<Type>::value &&		// Mat | Row | Col
		!std::is_same<typename Type::elem_type,qif::rat>::value
	>
> {

private:
	std::shared_ptr<Type> value_ref;

public:
	typedef typename Type::elem_type eT;
	static constexpr auto name = [] {
		if constexpr (std::is_same<Type, arma::Row<double>>::value)
			return _("Row<double>");
		else if constexpr (std::is_same<Type, arma::Col<double>>::value)
			return _("Col<double>");
		else if constexpr (std::is_same<Type, arma::Mat<double>>::value)
			return _("Mat<double>");
		else if constexpr (std::is_same<Type, arma::Row<uint>>::value)
			return _("Row<uint>");
		else if constexpr (std::is_same<Type, arma::Col<uint>>::value)
			return _("Col<uint>");
		else
			return _("Mat<uint>");
	}();

	bool load(handle src, bool convert) {
		const uint ndim = arma::is_Mat_only<Type>::value ? 2 : 1;

		using Array = array_t<eT, array::forcecast | array::f_style>;
		Array array;

		bool no_copy = isinstance<Array>(src);
		if(no_copy) {
			array = reinterpret_borrow<Array>(src);
			if(array.ndim() != ndim)
				return false;

			// We can avoid copy only if the data are contiguous and (for mat) in column-major order.
			// This can fail for matrices in row-major order, but also for vectors if they are created by slicing matrices.
			const ssize_t sz = sizeof(eT);
			no_copy = array.strides(0) == sz;

			if constexpr (arma::is_Mat_only<Type>::value) {
				no_copy = no_copy && array.strides(1) == sz * array.shape(0);
			}
		}

		// std::cout << "load double, covert " << convert << ", no_copy: " << no_copy << "\n";
		if(!no_copy) {
			// We need to copy/convert, fail if this is the no-convert pass
			if (!convert) return false;

			array = Array::ensure(src);
			if (!array || array.ndim() != ndim)
				return false;
			// std::cout << "array->mat with copy  " << array.data() << "\n";

			loader_life_support::add_patient(array);
		}

		if constexpr (arma::is_Mat_only<Type>::value) {
			value_ref.reset(new Type(
				(eT*)array.data(),
				array.shape(0),
				array.shape(1),
				false,  // access the same underlying memory
				true
			));
		} else {
			value_ref.reset(new Type(
				(eT*)array.data(),
				array.shape(0),
				false,  // access the same underlying memory
				true
			));
		}

		return true;
	}

    operator Type*() { return value_ref.get(); }
    operator Type&() { return *value_ref; }
    operator Type&&() && { return std::move(*value_ref); }
    template <typename _T> using cast_op_type = pybind11::detail::cast_op_type<_T>;

	// Normal returned non-reference, non-const value:
	static handle cast(Type&& src, return_value_policy /* policy */, handle /* parent */) {
		size_t sz = sizeof(eT);
		std::vector<size_t> shape, strides;
		// std::cout << "casting arma " << src.memptr() << "\n";

		if constexpr (arma::is_Mat_only<Type>::value) {
			shape = { src.n_rows, src.n_cols };
			strides = { sz, sz * src.n_rows };
		} else {
			shape = { src.n_elem };
			strides = { sz };
		}

		// move the array contents to a heap-based matrix, and create a capsule storing the pointer
		// while the array is being used
		Type *p = new Type(std::move(src));
		capsule base(p, [](void *o) {
			// std::cout << "deleting mat " << static_cast<Type*>(o)->memptr() << "\n";
			delete static_cast<Type *>(o);
		});

		return array(pybind11::dtype::of<eT>(), shape, strides, p->memptr(), base).release();
	}

    // const lvalue reference return; copy
	static handle cast(const Type& src, return_value_policy policy, handle parent) {
		return cast(Type(src), policy, parent);
	}
};


///////////////////////////// Mat<rat> | Row<rat> | Col<rat> ///////////////////////////////////////

template <typename Type>
class type_caster<
	Type,
	typename std::enable_if_t<
		arma::is_Mat<Type>::value &&		// Mat | Row | Col
		std::is_same<typename Type::elem_type,qif::rat>::value
	>
> {
	using rat = qif::rat;

	PYBIND11_TYPE_CASTER(Type, []{
		if constexpr (std::is_same<Type, arma::Row<rat>>::value)
			return _("Row<rat>");
		else if constexpr (std::is_same<Type, arma::Col<rat>>::value)
			return _("Col<rat>");
		else
			return _("Mat<rat>");
	}());

public:

	bool load(handle src, bool) {

		if(!isinstance<buffer>(src))
			return false;

		buffer buf = reinterpret_borrow<buffer>(src);
        buffer_info info = buf.request();
		if(info.format != "O")
			return false;		// dtype should be object


		if constexpr (std::is_same<Type, arma::Row<rat>>::value || std::is_same<Type, arma::Col<rat>>::value) {

			if(info.ndim != 1)
				return false;

			uint n_elem = info.shape[0];
			uint sx = info.strides[0];
			char* ptr = (char*)info.ptr;

			value.set_size(n_elem);

			for(uint x = 0; x < n_elem; x++)
				value(x) = handle(*(PyObject**)(ptr + x*sx)).cast<rat>();

		} else if constexpr (std::is_same<Type, arma::Mat<rat>>::value) {

			if(info.ndim != 2)
				return false;

			uint n_rows = info.shape[0];
			uint n_cols = info.shape[1];
			uint sx = info.strides[0];
			uint sy = info.strides[1];
			char* ptr = (char*)info.ptr;

			value.set_size(n_rows, n_cols);

			for(uint x = 0; x < n_rows; x++)
				for(uint y = 0; y < n_cols; y++)
					value(x,y) = handle(*(PyObject**)(ptr + x*sx + y*sy)).cast<rat>();
		}

		return true;
	}

    // const lvalue reference return; we always copy so no need to define other cases
	static handle cast(const Type& src, return_value_policy /* policy */, handle /* parent */) {

		size_t sz = sizeof(PyObject*);
		std::vector<size_t> shape, strides;

		if constexpr (arma::is_Mat_only<Type>::value) {
			shape = { src.n_rows, src.n_cols };
			strides = { sz, sz * src.n_rows };
		} else {
			shape = { src.n_elem };
			strides = { sz };
		}

		auto objs = new std::vector<PyObject*>();
		objs->reserve(src.n_elem);

		for(auto& q : src)
			objs->push_back(pybind11::cast(q).release().ptr());

		capsule base(objs, [](void *o) {
			auto objs = static_cast<std::vector<PyObject*>*>(o);
			for(auto obj : *objs)
				handle(obj).dec_ref();
			delete objs;
		});

		return array("O", shape, strides, &(*objs)[0], base).release();
	}
};

} // pybind11:detail