
namespace plot {

template <typename eT> using ME = lp::MatrixEntry<eT>;


template<typename eT>
void gnuplot_barycentric_3d(std::function<eT(Prob<eT>&)> f, std::string filename) {

	uint stepno = 30;

	std::ofstream file;
	file.open(filename);
	if(!file)
		throw std::runtime_error("cannot open " + filename);

	file << R"(
		# data in matrix form. The value for prior (a,b,c) should be in row a and column b
		set datafile separator ","
	)"
	<< "stepno = " << stepno << "\n"
	<< "$data << EOD\n";

	// The matrix sent to gnuplot has to cover the whole unit square. We could use NaN for the part
	// where a+b > 1, but this would cause the plot to be "cropped". Instead we limit b by setting b=1-a
	//
	for(uint i = 0; i <= stepno; i++) {
		for(uint j = 0; j <= stepno; j++) {
			uint jj = std::min(j, stepno - i);
			Prob<eT> pi = { eT(i), eT(jj), eT(stepno - i - jj) };
			pi /= stepno;
			file << f(pi) << (j < stepno ? "," : "\n");
		}
	}

	file << R"(EOD

		# compute the maximum function value in STATS_max. Needs to be done
		# before setting xrange/yrange cause stats filters the data
		stats $data matrix using (valid($3) ? $3 : 0) nooutput

		set terminal qt size 600,600 font 'Verdana,10'
		# set terminal pdf

		# use barycentric-coordinates, each prior (a,b,c) is displayed as a point in an
		# equilateral triangle. The cartesian coordinates (x,y) of the point are (c/2 + b, c * sqrt(3)/2)
		# (this choice makes the origin (0,0) coincite with the point prior on x1 (1,0,0).
		bary_x(a,b) = compl(a+b)/2 + b
		bary_y(a,b) = compl(a+b) * sqrt(3)/2
		compl(x) = 1.0 - x

		# To make the plot prettier we hide the cartesian axes and show the equilateral triangle and z-axes
		# on each of the 3 vertices.
		#
		unset border
		unset xtics
		unset ytics

		# show the triangle with vertices x_1, x_2, x_3
		set arrow 1 from bary_x(1,0),bary_y(1,0) to bary_x(0,0),bary_y(0,0) nohead front lt -1 lw 1
		set arrow 2 from bary_x(0,0),bary_y(0,0) to bary_x(0,1),bary_y(0,1) nohead front lt -1 lw 1
		set arrow 3 from bary_x(0,1),bary_y(0,1) to bary_x(1,0),bary_y(1,0) nohead front lt -1 lw 1

		set label 1 "x_1" at bary_x(1,0)+.02,bary_y(1,0)-.05		# small offset wrt the actual vertex
		set label 2 "x_2" at bary_x(0,1)+.02,bary_y(0,1)
		set label 3 "x_3" at bary_x(0,0),    bary_y(0,0)+0.02

		# show z-axis on each vertex. For x_1 we use the zzeroaxis that includes ticks. For x_2,x_3 we draw ourselves
		set zzeroaxis
		set arrow 4 from bary_x(0,1),bary_y(0,1),0 to bary_x(0,1),bary_y(0,1),STATS_max nohead front dt 3 lc rgb '#000000' lw 1
		set arrow 5 from bary_x(0,0),bary_y(0,0),0 to bary_x(0,0),bary_y(0,0),STATS_max nohead front dt 3 lc rgb '#000000' lw 1

		# ranges
		set xrange [0:1]
		set yrange [0:1]
		set zrange [0:]
		set ticslevel 0			# set the XY-plane at z=0. Redundant if zrange starts at 0

		# these affect how the plot is shown
		set view 61, 41			# rotate the plot
		set hidden3d
		# set pm3d
		# set palette gray
		# set contour
		# set cntrparam levels incremental -1, 0.2, 1

		# coordinate transformation: normalize with stepno, limit b if a+b > 1
		min(x,y) = x < y ? x : y
		xx(a,b) = bary_x(a/stepno, min(b, stepno-a)/stepno)
		yy(a,b) = bary_y(a/stepno, min(b, stepno-a)/stepno)

		# the actual plot. gnuplot uses $1 for the _column_ index of the matrix
		# and $2 for the row index, which is weird so we swap them
		splot $data matrix \
			using (xx($2,$1)):(yy($2,$1)):($3) \
			title '' with lines linestyle 1

		pause mouse close

	)";
	file.close();
}

} // namespace plot
