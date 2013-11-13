#include "LinearProgram.h"
#include <stdio.h>
#include <stdlib.h>
/*
This file belongs to the LIBQIF library.
A Quantitative Information Flow C++ Toolkit Library.
Copyright (C) 2013  Universidad Nacional de Río Cuarto(National University of Río Cuarto).
Author: Martinelli Fernán - fmartinelli89@gmail.com - Universidad Nacional de Río Cuarto (Argentina)
LIBQIF Version: 1.0
Date: 12th Nov 2013 
========================================================================
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

=========================================================================
*/
VectorType LinearProgram::solve(StringType& equality,StringType& inequality,StringType& objective){
	MatrixType eq=arma::mat(equality); // the matrix with the equalities
	MatrixType ineq=arma::mat(inequality); // the matrix with the inequalities
	VectorType obj=arma::vec(objective); //the objective function as a vector
	return solve(eq,ineq,obj);
}

VectorType LinearProgram::solve(StringType& equality,StringType& inequality,StringType& objective,StringType& rows_constraints){
	MatrixType eq=arma::mat(equality); // the matrix with the equalities
	MatrixType ineq=arma::mat(inequality); // the matrix with the inequalities
	VectorType obj=arma::vec(objective); //the objective function as a vector
	MatrixType constraints=arma::mat(rows_constraints); // the matrix with the lineal system constrains
	return solve(eq,ineq,obj,constraints);
}
	
VectorType LinearProgram::solve(MatrixType equality,MatrixType inequality,VectorType objective){
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
	for(i=1;i<=n_variables;++i){
		glp_set_col_bnds(lp, i, GLP_LO, 0.0, -inequality.at(i-1,inequality.n_cols-1));
	}
	//---------------------------------------------
	
	//objective coefficients-----------------------
	for(i=1;i<=n_variables;++i){
		glp_set_obj_coef(lp, i, objective.at(i-1));
	}
	//---------------------------------------------
	
	//Row indices of each element are stored in the array ia, column indices are stored in the array ja, and numerical
	//values of corresponding elements are stored in the array ar.
	int j;
	for(i=0;i<n_variables;++i){
		for(j=0;j<n_constraints;++j){
			ia[i*n_constraints+j]=i+1, ja[i*n_constraints+j]=j+1, ar[i*n_constraints+j]=equality.at(i,j);
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
	for(i=0;i<n_variables;++i){
		x[i] = glp_get_col_prim(lp, i+1);
	}
	
	printf("\nz = %g; x1 = %g; x2 = %g; x3 = %g\n",z, x[1], x[2], x[3]);


	//frees all the memory allocated to the problem object
	glp_delete_prob(lp);

	//Its not implemented yet
	VectorType a;
	return a;
}

VectorType LinearProgram::solve(MatrixType equality,MatrixType inequality,VectorType objective,MatrixType rows_constraints){
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
	for(i=0;i<rows_constraints.n_rows;++i){
		glp_set_row_bnds(lp, rows_constraints.at(i,0), GLP_UP, 0.0, -inequality.at(i,inequality.n_cols-1));
	}
	//---------------------------------------------
	
	//cols: x1..xn---------------------------------
	glp_add_cols(lp, n_variables);
	
	//bounds : ineq
	for(i=1;i<=n_variables;++i){
		glp_set_col_bnds(lp, i, GLP_LO, 0.0, -inequality.at(i-1,inequality.n_cols-1));
	}
	//---------------------------------------------
	
	//objective coefficients-----------------------
	for(i=1;i<=n_variables;++i){
		glp_set_obj_coef(lp, i, objective.at(i-1));
	}
	//---------------------------------------------
	
	//Row indices of each element are stored in the array ia, column indices are stored in the array ja, and numerical
	//values of corresponding elements are stored in the array ar.
	int j;
	for(i=0;i<n_variables;++i){
		for(j=0;j<n_constraints;++j){
			ia[i*n_constraints+j]=i+1, ja[i*n_constraints+j]=j+1, ar[i*n_constraints+j]=equality.at(i,j);
		}
	}
	//calls the routine glp_load_matrix, which loads information from these three arrays into the problem object.
	glp_load_matrix(lp, n_variables*n_constraints, ia, ja, ar);
	//calls the routine glp_simplex, which is a driver to the simplex method, in order to solve the LP problem. 
	//This routine finds an optimal solution and stores all relevant information back into the problem object.
	glp_simplex(lp, NULL);
	//obtains a computed value of the objective function
	z=glp_get_obj_val(lp);
	//obtain computed values of structural variables (columns), 
	//which correspond to the optimal basic solution found by the solver.
	for(i=0;i<n_variables;++i){
		x[i]=glp_get_col_prim(lp, i+1);
	}
	
	printf("\nz = %g; x1 = %g; x2 = %g; x3 = %g\n",z, x[1], x[2], x[3]);


	//frees all the memory allocated to the problem object
	glp_delete_prob(lp);

	//Its not implemented yet
	VectorType a;
	return a;
}

