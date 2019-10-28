
namespace channel::comp {

template<typename eT>
inline
Chan<eT> parallel(const Chan<eT>& C1, const Chan<eT>& C2) {
	if(C1.n_rows != C2.n_rows)
		throw std::runtime_error("rows mismatch");

	Chan<eT> C(C1.n_rows, C1.n_cols * C2.n_cols);
	for(uint i = 0; i < C1.n_rows; i++)
	for(uint j = 0; j < C1.n_cols; j++)
	for(uint k = 0; k < C2.n_cols; k++)
		C(i, j*C2.n_cols + k) = C1(i,j) * C2(i,k);

	return C;
}

template<typename eT>
inline
Chan<eT> repeated_independent(const Chan<eT>& C, uint n) {
	if(n == 0)
		return channel::no_interference<eT>(C.n_rows);
	if(n == 1)
		return C;

	Chan<eT> Cn = parallel(C, C);
	for(n -= 2; n > 0; n--)
		Cn = parallel(Cn, C);

	return Cn;
}


} // namespace channel::comp
