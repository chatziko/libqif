
namespace games {

using std::vector;


template<typename eT>
eT
bayes_subgradient(const Prob<eT>& pi, const vector<vector<Chan<eT>>>& Cs, const Prob<eT>& delta, arma::Row<eT>& subgrad) {
	uint n_adv = Cs.size();
	uint n_def = delta.n_elem;
	uint n_rows = Cs[0][0].n_rows;
	uint n_cols = Cs[0][0].n_cols;

	eT max(-1);

	for(uint a = 0; a < n_adv; a++) {	// max_a

		arma::Row<eT> subgrad_temp(n_def);
		subgrad_temp.fill(0);
		eT sumy(0);

		for(uint y = 0; y < n_cols; y++) { //sum_y

			arma::Row<eT> subgrad2;
			eT max2(-1);

			for(uint x = 0; x < n_rows; x++) {
				arma::Row<eT> subgrad2_temp(n_def);
				for(uint d = 0; d < n_def; d++)
					subgrad2_temp(d) = pi(x) * Cs[a][d](x, y);

				eT z = arma::cdot(delta, subgrad2_temp);
				if(z > max2) {
					max2 = z;
					subgrad2 = subgrad2_temp;
				}
			} // max_x

			subgrad_temp += subgrad2;
			sumy += max2;
		} // sum_y

		if(sumy > max) {
			subgrad = subgrad_temp;
			max = sumy;
		}
	} // max_a

	return max;
}

template<typename eT>
std::pair<eT, Prob<eT>>
minmax_hidden_bayes(const Prob<eT>& pi, const vector<vector<Chan<eT>>>& Cs) {

	uint n_def = Cs[0].size();

	// Notes secton 3.1: R should be such that R >= ||x^1 - x^*||_2. We start from a
	// uniform, and the farthest we can go is a point distribution, in which case
	// R = sqrt( n-1 / n^2 + 1 - 1/n ) = sqrt( 1 - 1/n )
	//
	const eT R2 = eT(1) - 1/eT(n_def);		// R^2

	Prob<eT> delta = probab::uniform<eT>(n_def);		// start from uniform
	Prob<eT> delta_best;
	eT f_best = infinity<eT>();
	eT l_best(0);
	eT sum1(0), sum2(0), sum3(0);

	for(eT k = 1; true; k++) {
//		eT alpha = eT(1) / k;						// square summable, see Boyd page 3
		eT alpha = eT(0.1) / sqrt(k);				// Nonsummable diminishing, see Boyd page 3
		arma::Row<eT> g;
		eT f = bayes_subgradient(pi, Cs, delta, g);	// returns f, stores subgrad in g

		// keep the best f/delta
		if(f < f_best) {
			f_best = f;
			delta_best = delta;
		}

		// stopping criterion, see notes Section 3.4
		sum1 += alpha * f;
		sum2 += alpha * alpha * arma::cdot(g, g);
		sum3 += alpha;
		eT l = (2 * sum1 - R2 - sum2) / (2 * sum3);
		if(l > l_best) {
			l_best = l;
			std::cerr << "\r" << f_best << "   " << (f_best - l_best) << "                 ";
		}

		if(f_best - l_best < 1e-2)
			break;

		// update and project back to the probability simplex
		delta -= alpha * g;
		delta = metric::optimize::simplex_project(delta);
	}
	std::cerr << "\r";

	return std::pair<eT, Prob<eT>>(f_best, delta_best);
}



template<typename eT>
eT
bayes_both_subgradients(
	const Prob<eT>& pi,
	const vector<vector<Chan<eT>>>& Cs,
	const Prob<eT>& alpha,
	const Prob<eT>& delta,
	arma::Row<eT>& subgrad_a,
	arma::Row<eT>& subgrad_d
) {
	uint n_adv = alpha.n_elem;
	uint n_def = delta.n_elem;
	uint n_rows = Cs[0][0].n_rows;
	uint n_cols = Cs[0][0].n_cols;

	subgrad_a.set_size(n_adv);
	subgrad_d.set_size(n_def);
	subgrad_a.fill(0);
	subgrad_d.fill(0);

	eT f(0);

	for(uint a = 0; a < n_adv; a++) {	// sum_a

		for(uint y = 0; y < n_cols; y++) { //sum_y

			arma::Row<eT> max_sg;
			eT max(-1);
			arma::Row<eT> max_sg_temp(n_def);

			for(uint x = 0; x < n_rows; x++) {
				for(uint d = 0; d < n_def; d++)
					max_sg_temp(d) = pi(x) * Cs[a][d](x, y);

				eT z = arma::cdot(delta, max_sg_temp);
				if(z > max) {
					max = z;
					max_sg = max_sg_temp;
				}
			} // max_x

			subgrad_a(a) += max;
			subgrad_d += alpha(a) * max_sg;
			f += alpha(a) * max;

		} // sum_y

	} // sum_a

	return f;
}

template<typename eT>
std::pair<eT, Prob<eT>>
minmax_hidden_bayes_both(
	const Prob<eT>& pi,
	const vector<vector<Chan<eT>>>& Cs
) {
	uint n_adv = Cs.size();
	uint n_def = Cs[0].size();

	// see minmax_hidden_bayes
	const eT R2 = eT(1) - 1/eT(n_def);		// R^2

	Prob<eT> alpha = probab::uniform<eT>(n_adv);		// start from uniform
	Prob<eT> delta = probab::uniform<eT>(n_def);		// start from uniform
	Prob<eT> alpha_avg(n_adv);
	Prob<eT> delta_avg(n_def);
	LargeAvg<eT> f_avg;
	eT l_best(0);
	vector<LargeAvg<eT>> alpha_lavg(n_adv);
	vector<LargeAvg<eT>> delta_lavg(n_adv);
	eT sum1(0), sum2(0), sum3(0);

	for(eT k = 1; true; k++) {
//		eT gamma = eT(1) / k;						// square summable, see Boyd page 3
		eT gamma = eT(0.1) / sqrt(k);				// Nonsummable diminishing, see Boyd page 3
		// eT gamma = eT(1);
		arma::Row<eT> g_a, g_d;

		eT f = bayes_both_subgradients(pi, Cs, alpha, delta, g_a, g_d);	// stores subgrads in g_a/g_d

		// keep track over the average f/alpha/delta
		f_avg.add(f);

		for(uint a = 0; a < n_adv; a++)
			alpha_avg(a) = alpha_lavg[a].add(alpha(a));

		for(uint d = 0; d < n_def; d++)
			delta_avg(d) = delta_lavg[d].add(delta(d));

		// print (both f_avg and f(alpha_avg, delta_avg)
		arma::Row<eT> g1, g2;		// unused, just to call bayes_both_subgradients
		eT ff = bayes_both_subgradients(pi, Cs, alpha_avg, delta_avg, g1, g2);
		std::cout << "\r" << f_avg.value() << ", " << ff << ", " << std::abs(f_avg.value() - 0.087649797994) << "     ";

		// stopping criterion, see notes Section 3.4
		// sum1 += gamma * f;
		// sum2 += gamma * gamma * arma::cdot(g_d, g_d);
		// sum3 += gamma;
		// eT l = (2 * sum1 - R2 - sum2) / (2 * sum3);
		// if(l > l_best) {
		// 	l_best = l;
		// 	std::cerr << "\r" << f_avg.value() << "   " << (f_avg.value() - l_best) << "                 ";
		// }

		// if(f_avg.value() - l_best < 1e-2)
		// 	break;

		// update and project back to the probability simplex
		alpha += gamma * g_a;
		alpha = metric::optimize::simplex_project(alpha);

		delta -= gamma * g_d;
		delta = metric::optimize::simplex_project(delta);
	}
	std::cerr << "\r";

	return std::pair<eT, Prob<eT>>(f_avg.value(), delta_avg);	// TODO return also alpha_avg?
}

template<typename eT>
std::pair<eT, Prob<eT>>
minmax_hidden_bayes_both2(
	const Prob<eT>& pi,
	const vector<vector<Chan<eT>>>& Cs
) {
	uint n_adv = Cs.size();
	uint n_def = Cs[0].size();

	Prob<eT> alpha = probab::uniform<eT>(n_adv);		// start from uniform
	Prob<eT> delta = probab::uniform<eT>(n_def);		// start from uniform
	Prob<eT> delta1 = probab::uniform<eT>(n_def);		// start from uniform
	Prob<eT> alpha_avg(n_adv);
	Prob<eT> delta_avg(n_def);
	LargeAvg<eT> f_avg;
	eT f_best = infinity<eT>();
	vector<LargeAvg<eT>> alpha_lavg(n_adv);
	vector<LargeAvg<eT>> delta_lavg(n_adv);

	for(eT k = 1; true; k++) {
//		eT gamma = eT(1) / k;						// square summable, see Boyd page 3
		eT gamma = eT(0.1) / sqrt(k);				// Nonsummable diminishing, see Boyd page 3
		// eT gamma = eT(1);
		arma::Row<eT> g;
		arma::Row<eT> g_a, g_d;

		// first our method
		eT f = bayes_subgradient(pi, Cs, delta1, g);	// returns f, stores subgrad in g
		if(f < f_best) {
			f_best = f;
		}

		// then theirs
		f = bayes_both_subgradients(pi, Cs, alpha, delta, g_a, g_d);	// stores subgrads in g_a/g_d
		f_avg.add(f);

		std::cerr << "\r" << f_best << ", " << f_avg.value() << "   " << std::abs(f_avg.value() - f_best) << "                 ";

		// update 1
		delta1 -= gamma * g;
		delta1 = metric::optimize::simplex_project(delta1);

		// update and project back to the probability simplex
		alpha += gamma * g_a;
		alpha = metric::optimize::simplex_project(alpha);

		delta -= gamma * g_d;
		delta = metric::optimize::simplex_project(delta);
	}
	std::cerr << "\r";

	return std::pair<eT, Prob<eT>>(f_avg.value(), delta_avg);	// TODO return also alpha_avg?
}


template<typename eT>
std::pair<eT, Prob<eT>>
dp_hidden_lp(const vector<Chan<eT>>& Cs) {

	uint n_def = Cs.size();
	uint n_rows = Cs[0].n_rows;
	uint n_cols = Cs[0].n_cols;

	// start with uniform dist
	Prob<eT> delta = probab::uniform<eT>(n_def);
	Chan<eT> C(n_rows, n_cols);

	while(true) {
		// create current matrix
		C.fill(eT(0));
		for(uint i = 0; i < n_def; i++)
			C += delta(i) * Cs[i];

		// compute lambda
		eT lambda(0);
		for(uint x1 = 0; x1 < n_rows; x1++)
		for(uint x2 = 0; x2 < n_rows; x2++)
		for(uint y = 0; y < n_cols; y++)
			lambda = max(lambda, C(x1,y) / C(x2,y));

		// solve lp	
		lp::LinearProgram<eT> lp;

		auto z = lp.make_var(-infinity<eT>());
		auto delta_v = lp.make_vars(n_def, eT(0));

		// minimize z
		lp.maximize = false;
		lp.set_obj_coeff(z, 1);

		// sum_i delta(i) = 1
		auto con = lp.make_con(eT(1), eT(1));
		for(uint i = 0; i < n_def; i++)
			lp.set_con_coeff(con, delta_v[i], eT(1));

		// for each x,x',y
		//   z >= sum_i C_i(x,y) gamma(i)  - lambda( sum_i C_i(x',y) gamma(i) )      <=>
		//   z >= sum_i [C_i(x,y) - lambda C_i(x',y)] gamma(i)
		//
		for(uint x1 = 0; x1 < n_rows; x1++) {
		for(uint x2 = 0; x2 < n_rows; x2++) {
			if(x1 == x2)
				continue;

			for(uint y = 0; y < n_cols; y++) {
				auto con = lp.make_con(-infinity<eT>(), eT(0));

				lp.set_con_coeff(con, z, eT(-1));

				for(uint i = 0; i < n_def; i++)
					lp.set_con_coeff(con, delta_v[i], Cs[i](x1,y) - lambda * Cs[i](x2, y));
			}
		}}

		// ready to solve
		lp.msg_level = lp::MsgLevel::ON;
		if(!lp.solve())
			throw std::runtime_error("problem should be solvable");
		std::cout << "SOLVED: lambda: " << lambda << ", F(lambda):" << lp.objective() << "\n";

		for(uint i = 0; i < n_def; i++)
			delta(i) = lp.solution(delta_v[i]);

		// stop if the objective function reaches 0
		if(equal(lp.objective(), eT(0), eT(1e-5)))
			return std::pair<eT, Prob<eT>>(lambda, delta);
	}
}

template<typename eT>
std::tuple<eT, uint, Prob<eT>>
dp_hidden(const vector<vector<Chan<eT>>>& Cs) {

	auto [ max, delta_max ] = dp_hidden_lp(Cs[0]);
	uint alpha_max = 0;

	uint n_adv = Cs.size();
	for(uint a = 1; a < n_adv; a++) {
		auto [ res, delta ] = dp_hidden_lp(Cs[a]);

		if(less_than(max, res)) {
			max = res;
			alpha_max = a;
			delta_max = delta;
		}
	}

	return std::tuple<eT, uint, Prob<eT>>(max, alpha_max, delta_max);
}




} // namespace games
