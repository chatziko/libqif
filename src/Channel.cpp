#include "Channel.h"
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
Channel::Channel (double** new_channel )
{
	//VER COMO IMPLEMENTAR ESTO.. NO ESTA DIRECTAMENTE EN ARMADILLO
}
*/
Channel::Channel(StringType& new_channel_elements)
{
	matrix=arma::mat(new_channel_elements);
	str=new_channel_elements;
	if(!this->rep_ok()){throw 1;}
}

Channel::Channel(MatrixType& new_channel_elements)
{
	matrix=arma::mat(new_channel_elements);
	if(!this->rep_ok()){throw 1;}
}
/*
Channel::~Channel()
{
	matrix.~mat();
}
*/

Channel Channel::new_id_channel(IntType size)
{
	if(size<0){throw 1;}
	std::string new_channel_elements=std::string("");
	int i,j;
	for(i=0;i<size;++i){
		for(j=0;j<size;++j){
			if(i==j)
				new_channel_elements+=std::string(" 1");
			else
				new_channel_elements+=std::string(" 0");
		}
		new_channel_elements+=std::string(" ; ");	
	}
	str=new_channel_elements;
	return Channel(new_channel_elements);
}

Channel Channel::clone()
{
	return Channel(matrix);
}

int Channel::inputs_number()
{
	return matrix.n_rows;
}

int Channel::outputs_number()
{
	return matrix.n_cols;
}

bool Channel::is_symmetric()
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
bool Channel::is_partial_symmetric()
{
	//Its not implemented yet
	return true;
}

bool Channel::is_equal_to(const Channel& other)
{
	//Its not implemented yet
	return true;
}
*/
VectorType Channel::get_row(IntType index)
{
	return matrix.row(index);
}
		
VectorType Channel::get_column(IntType index)
{
	return matrix.col(index);
}
/*		
void Channel::set_row(IntType index,VectorType& new_row_elements)
{
	
}

void Channel::set_row(IntType index,StringType& new_row_elements)
{
	arma::vec new_row= arma::vec(new_row_elements);
	set_row(index,new_row);
}
*/
DoubleType Channel::at(IntType index_x,IntType index_y)
{
	return matrix.at(index_x,index_y);
}

