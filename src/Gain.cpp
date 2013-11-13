#include "Gain.h"
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
/*
Gain::Gain (double** new_gain_function )
{
	//Its not implemented yet
	//int new_size= sizeof(new_gain_function)/sizeof(*new_gain_function);
}
*/
Gain::Gain(StringType& new_gain_elements)
{
	matrix=arma::mat(new_gain_elements);
	str=new_gain_elements;
	//if(!this->rep_ok()){throw 1;}
}

Gain::Gain(MatrixType& new_gain_elements)
{
	matrix=arma::mat(new_gain_elements);
	//if(!this->rep_ok()){throw 1;}
}
/*
Gain::~Gain()
{
	matrix.~mat();
}
*/
int Gain::inputs_number()
{
	return matrix.n_cols;
}

int Gain::guesses_number()
{
	return matrix.n_rows;
}
/*
Gain Gain::zeros(IntType inputs,IntType outputs)
{
	matrix= arma::mat(inputs,outputs);
	matrix.zeros();
	return Gain(matrix);
}

Gain Gain::ones(IntType inputs,IntType outputs)
{
	matrix=arma::mat(inputs,outputs);
	matrix.ones();
	return Gain(matrix);
}
*/
Gain Gain::new_id_function(IntType size)
{
	if(size<0){throw 1;}
	std::string new_gain_elements=std::string("");
	int i,j;
	for(i=0;i<size;++i){
		for(j=0;j<size;++j){
			if(i==j)
				new_gain_elements+=std::string(" 1");
			else
				new_gain_elements+=std::string(" 0");
		}
		new_gain_elements+=std::string(" ; ");	
	}
	str=new_gain_elements;
	return Gain(new_gain_elements);
}

Gain Gain::clone()
{
	return Gain(matrix);
}

bool Gain::is_symmetric()
{
	bool flag=(matrix.n_cols == matrix.n_rows);
	int i,j;
	for(i=0;i<matrix.n_cols && flag;++i){
		for(j=0;j<matrix.n_rows && flag;++j){
			flag=flag && matrix.at(i,j)==matrix.at(j,i);
		}
	}
	return flag;
}
/*
bool Gain::is_partial_symmetric()
{
	//Its not implemented yet
	return true;
}

bool Gain::is_equal_to(const Gain& other)
{
	bool flag= (matrix.n_rows == other.n_rows);
	flag= flag && (matrix.n_cols == other.n_cols);
	int i,j;
	for(i=0;i<matrix.n_cols && flag;i++){
		for(j=0;j<matrix.n_rows && flag;j++){
			flag= flag && matrix.at(i,j)==other.at(j,i);
		}
	}
	return flag;
	return true;
}*/

VectorType Gain::get_row(IntType index)
{
	return matrix.row(index);
}
		
VectorType Gain::get_column(IntType index)
{
	return matrix.col(index);
}
/*		
void Gain::set_row(IntType index,arma::vec& new_row_elements)
{
	
}

void Gain::set_row(IntType index,std::string& new_row_elements)
{
	arma::vec new_row= arma::vec(new_row_elements);
	set_row(index,new_row);
}
*/
DoubleType Gain::at(IntType index_x ,IntType index_y){
	return matrix.at(index_x,index_y);
}