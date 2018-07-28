// EnableIf alias for SFINAE to make things more readable.
// We use the one but last method from http://loungecpp.wikidot.com/tips-and-tricks:enable-if-for-c-11
// This allows to have functions with exactly the same signature, differing only in the EnableIf/DisableIf conditions
// Note: the very last method with parameter packs (i.e. EnableIf<cond>...) does not work with clang (https://llvm.org/bugs/show_bug.cgi?id=11723)
//
enum class enabled { _ };
constexpr auto _ = enabled::_;
template <typename T> using Invoke = typename T::type;
template <typename Condition> using EnableIf  = Invoke<std::enable_if<  Condition::value, enabled >>;
template <typename Condition> using DisableIf = Invoke<std::enable_if< !Condition::value, enabled >>;

typedef uint32_t uint;
typedef mpq_class rat;

// default eT type, allow overriding before including qif
#ifndef QIF_DEFAULT_ET
#define QIF_DEFAULT_ET double
#endif
typedef QIF_DEFAULT_ET eT_def;

// Chan, Prob
using arma::Mat;
using arma::mat;
using arma::Row;
using arma::Col;

template<typename eT> using Chan = Mat<eT>;
template<typename eT> using Prob = Row<eT>;

typedef Chan<double> chan;
typedef Chan<float> fchan;
typedef Chan<rat>   rchan;

typedef Row<double>  prob;
typedef Row<float>  fprob;
typedef Row<rat>    rprob;

typedef Mat<rat>     rmat;
typedef Col<rat>  rcolvec;
typedef Row<rat>  rrowvec;

// Metric
// R: result type, what we measure distances in
// T: metric space type, what elements we measure the distance of
template<typename R, typename T>
class Metric : public std::function<R(const T&, const T&)> {
	public:
	std::function<bool(const T&, const T&)> chainable =
		[](const T&, const T&) -> bool { return false; };		// safe default, see metric.h

	using std::function<R(const T&, const T&)>::function;
};

// Point
template<typename eT>
struct Point {
	eT x, y;
	Point() {}
	Point(eT x, eT y) : x(x), y(y) {}

	bool operator==(const Point<eT>& rhs) const { return equal(this->x, rhs.x) && equal(this->y, rhs.y); }
	Point<eT> operator+(const Point<eT>& rhs) const { return Point<eT>(this->x + rhs.x, this->y + rhs.y); }
};

typedef Point<double> point;
typedef Point<float> fpoint;
typedef Point<rat>   rpoint;

// LatLog
template<typename eT>
struct LatLon {
	eT lat, lon;
	LatLon() {}
	LatLon(eT lat, eT lon) : lat(lat), lon(lon) {}

	LatLon<eT> add_vector(eT distance, eT angle);
};

typedef LatLon<double> latlon;
typedef LatLon<float> flatlon;
typedef LatLon<rat>   rlatlon;

// templates for SFINAE
template<typename T>  using is_Chan = arma::is_Mat_only<T>;
template<typename T>  using is_Prob = arma::is_Row<T>;

template<typename T>  struct is_Point            { static const bool value = false; };
template<typename eT> struct is_Point<Point<eT>> { typedef eT elem_type; static const bool value = true;  };
