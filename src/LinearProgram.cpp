#include "LinearProgram.h"
#include <stdio.h>
#include <stdlib.h>
/* Author: Fernan Martinelli*/
vec LinearProgram::solve(char* equality,char * inequality,char * objective){
	return solve(mat(equality),mat(inequality),vec(objective));
}

vec LinearProgram::solve(char* equality,char * inequality,char * objective,char * rows_constraints){
	return solve(mat(equality),mat(inequality),vec(objective),mat(rows_constraints));
}
	
vec LinearProgram::solve(mat equality,mat inequality,vec objective){
	int n_variables= objective.size();
	int n_constraints = equality.n_rows;
	
	glp_prob *lp;
	int ia[1+1000], ja[1+1000];
	double ar[1+1000];
	double x[n_variables];
	double z;

	//creates a problem object
	lp = glp_create_prob();
	//calls the routine glp_set_obj_dir in order to set the optimization direction flag, where GLP_MAX means maximization.
	glp_set_obj_dir(lp, GLP_MAX);
	
	int i;
	//rows: eq ------------------------------------
	glp_add_rows(lp, n_constraints);
	//for(i=1;i<=n_constraints;i++){
	//	glp_set_row_bnds(lp, i, GLP_UP, 0.0, -inequality.at(i-1,inequality.n_cols-1));
	//}
	//---------------------------------------------
	
	//cols: x1..xn---------------------------------
	glp_add_cols(lp, n_variables);
	
	//bounds : ineq
	for(i=1;i<=n_variables;i++){
		glp_set_col_bnds(lp, i, GLP_LO, 0.0, -inequality.at(i-1,inequality.n_cols-1));
	}
	//---------------------------------------------
	
	//objective coefficients-----------------------
	for(i=1;i<=n_variables;i++){
		glp_set_obj_coef(lp, i, objective.at(i-1));
	}
	//---------------------------------------------
	
	//Row indices of each element are stored in the array ia, column indices are stored in the array ja, and numerical
	//values of corresponding elements are stored in the array ar.
	int j;
	for(i=0;i<n_variables;i++){
		for(j=0;j<n_constraints;j++){
			ia[i*n_constraints+j] = i+1, ja[i*n_constraints+j] = j+1, ar[i*n_constraints+j] = equality.at(i,j);
		}
	}
	//calls the routine glp_load_matrix, which loads information from these three arrays into the problem object.
	glp_load_matrix(lp, n_variables*n_constraints, ia, ja, ar);
	//calls the routine glp_simplex, which is a driver to the simplex method, in order to solve the LP problem. 
	//This routine finds an optimal solution and stores all relevant information back into the problem object.
	glp_simplex(lp, NULL);
	//obtains a computed value of the objective function
	z = glp_get_obj_val(lp);
	//obtain computed values of structural variables (columns), 
	//which correspond to the optimal basic solution found by the solver.
	for(i=0;i<n_variables;i++){
		x[i] = glp_get_col_prim(lp, i+1);
	}
	
	printf("\nz = %g; x1 = %g; x2 = %g; x3 = %g\n",z, x[1], x[2], x[3]);


	//frees all the memory allocated to the problem object
	glp_delete_prob(lp);

	//Its not implemented yet
	vec a;
	return a;
}

vec LinearProgram::solve(mat equality,mat inequality,vec objective,mat rows_constraints){
	int n_variables= objective.size();
	int n_constraints = equality.n_rows;
	
	glp_prob *lp;
	int ia[1+1000], ja[1+1000];
	double ar[1+1000];
	double x[n_variables];
	double z;

	//creates a problem object
	lp = glp_create_prob();
	//calls the routine glp_set_obj_dir in order to set the optimization direction flag, where GLP_MAX means maximization.
	glp_set_obj_dir(lp, GLP_MAX);
	
	int i;
	//rows: eq ------------------------------------
	glp_add_rows(lp, n_constraints);
	for(i=0;i<rows_constraints.n_rows;i++){
		glp_set_row_bnds(lp, rows_constraints.at(i,0), GLP_UP, 0.0, -inequality.at(i,inequality.n_cols-1));
	}
	//---------------------------------------------
	
	//cols: x1..xn---------------------------------
	glp_add_cols(lp, n_variables);
	
	//bounds : ineq
	for(i=1;i<=n_variables;i++){
		glp_set_col_bnds(lp, i, GLP_LO, 0.0, -inequality.at(i-1,inequality.n_cols-1));
	}
	//---------------------------------------------
	
	//objective coefficients-----------------------
	for(i=1;i<=n_variables;i++){
		glp_set_obj_coef(lp, i, objective.at(i-1));
	}
	//---------------------------------------------
	
	//Row indices of each element are stored in the array ia, column indices are stored in the array ja, and numerical
	//values of corresponding elements are stored in the array ar.
	int j;
	for(i=0;i<n_variables;i++){
		for(j=0;j<n_constraints;j++){
			ia[i*n_constraints+j] = i+1, ja[i*n_constraints+j] = j+1, ar[i*n_constraints+j] = equality.at(i,j);
		}
	}
	//calls the routine glp_load_matrix, which loads information from these three arrays into the problem object.
	glp_load_matrix(lp, n_variables*n_constraints, ia, ja, ar);
	//calls the routine glp_simplex, which is a driver to the simplex method, in order to solve the LP problem. 
	//This routine finds an optimal solution and stores all relevant information back into the problem object.
	glp_simplex(lp, NULL);
	//obtains a computed value of the objective function
	z = glp_get_obj_val(lp);
	//obtain computed values of structural variables (columns), 
	//which correspond to the optimal basic solution found by the solver.
	for(i=0;i<n_variables;i++){
		x[i] = glp_get_col_prim(lp, i+1);
	}
	
	printf("\nz = %g; x1 = %g; x2 = %g; x3 = %g\n",z, x[1], x[2], x[3]);


	//frees all the memory allocated to the problem object
	glp_delete_prob(lp);

	//Its not implemented yet
	vec a;
	return a;
}

