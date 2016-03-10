using arma::Mat;
using arma::mat;
using arma::Row;
using arma::Col;

// helper aliases for SFINAE, see http://loungecpp.wikidot.com/tips-and-tricks:enable-if-for-c-11
enum class enabled {}; // just a type that can be used as a template parameter and is as inocuous as possible
template <typename T> using Invoke = typename T::type;
template <typename Condition> using EnableIf  = Invoke<std::enable_if< Condition::value, enabled>>;

template <typename T> using eT = typename T::elem_type;

template<typename eT> using Chan    = Mat<eT>;
template<typename T>  using is_Chan = arma::is_Mat_only<T>;

template<typename eT> using Prob    = Row<eT>;
template<typename T>  using is_Prob = arma::is_Row<T>;

typedef uint32_t uint;

typedef mpq_class rat;

typedef Chan<double> chan;
typedef Chan<float> fchan;
typedef Chan<rat>   rchan;

typedef Row<double>  prob;
typedef Row<float>  fprob;
typedef Row<rat>    rprob;

typedef Mat<rat>     rmat;
typedef Col<rat>  rcolvec;
typedef Row<rat>  rrowvec;


template<typename eT>
struct Point {
	eT x, y;
	Point() {}
	Point(eT x, eT y) : x(x), y(y) {}

	bool operator==(const Point<eT>& rhs) const { return equal(this->x, rhs.x) && equal(this->y, rhs.y); }
};

template<typename T>
struct is_Point { static const bool value = false; };
template<typename eT>
struct is_Point<Point<eT>> { typedef eT elem_type; static const bool value = true;  };

typedef Point<double> point;
typedef Point<float> fpoint;
typedef Point<rat>   rpoint;

