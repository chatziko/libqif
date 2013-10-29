#include "Channel.h"
/* Author: Fernan Martinelli*/

/*
Channel::Channel (double** new_channel )
{
	//VER COMO IMPLEMENTAR ESTO.. NO ESTA DIRECTAMENTE EN ARMADILLO
}
*/
Channel::Channel (std::string& new_channel_elements )
{
	C=mat(new_channel_elements);
}

Channel::Channel(mat new_channel_elements )
{
	C=mat(new_channel_elements);
}

Channel::~Channel()
{
	C.~mat();
}

Channel Channel::new_id_channel ( int size )
{
	std::string new_channel_elements=std::string("");
	int i,j;
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			if(i==j)
				new_channel_elements+=std::string(" 1");
			else
				new_channel_elements+=std::string(" 0");
		}
		new_channel_elements+=std::string(" ; ");	
	}
	return Channel(new_channel_elements);
}

Channel Channel::clone()
{
	return Channel(C);
}

int Channel::inputs_number()
{
	return C.n_rows;
}

int Channel::outputs_number()
{
	return C.n_cols;
}

bool Channel::is_symmetric()
{
	bool flag= (C.n_cols == C.n_rows);
	int i,j;
	for(i=0;i<C.n_cols && flag;i++){
		for(j=0;j<C.n_rows && flag;j++){
			flag= flag && C.at(i,j)==C.at(j,i);
		}
	}
	return flag;
}

bool Channel::is_partial_symmetric()
{
	//Its not implemented yet
	return true;
}

bool Channel::is_equal_to ( const Channel& other )
{
	//Its not implemented yet
	return true;
}

vec Channel::get_row ( int index )
{
	return C.row(index);
}
		
vec Channel::get_column ( int index ){
	return C.col(index);
}
		
void Channel::set_row ( int index,vec new_row_elements ){
	//Its not implemented yet
}

void Channel::set_row ( int index,char* new_row_elements ){
	vec new_row= vec(new_row_elements);
	set_row (index,new_row);
}

double Channel::at ( int index_x ,int index_y ){
	return C.at(index_x,index_y);
}
