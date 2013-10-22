#include "Gain.h"
#include <string>

/*
Gain::Gain (double** new_gain_function )
{
	//Its not implemented yet
	//int new_size= sizeof(new_gain_function)/sizeof(*new_gain_function);
}
*/
Gain::Gain (char* new_gain_elements )
{
	C=mat(new_gain_elements);
}

Gain::Gain (mat new_gain_elements )
{
	C=mat(new_gain_elements);
}

Gain::~Gain()
{
	C.~mat();
}

int Gain::inputs_number()
{
	return C.n_cols;
}

int Gain::guesses_number()
{
	return C.n_rows;
}

Gain Gain::zeros ( int inputs, int outputs )
{
	C= mat(inputs,outputs);
	C.zeros();
	return Gain(C);
}

Gain Gain::ones ( int inputs, int outputs )
{
	C= mat(inputs,outputs);
	C.ones();
	return Gain(C);
}

Gain Gain::new_id_function ( int size )
{
	std::string new_gain_elements=std::string("");
	int i,j;
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			if(i==j)
				new_gain_elements+=std::string(" 1");
			else
				new_gain_elements+=std::string(" 0");
		}
		new_gain_elements+=std::string(" ; ");	
	}
	return Gain(new_gain_elements.c_str());
}

Gain Gain::clone()
{
	return Gain(C);
}

bool Gain::is_symmetric()
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

bool Gain::is_partial_symmetric()
{
	//Its not implemented yet
	return true;
}

bool Gain::is_equal_to (const Gain& other )
{
	/*bool flag= (C.n_rows == other.n_rows);
	flag= flag && (C.n_cols == other.n_cols);
	int i,j;
	for(i=0;i<C.n_cols && flag;i++){
		for(j=0;j<C.n_rows && flag;j++){
			flag= flag && C.at(i,j)==other.at(j,i);
		}
	}
	return flag;*/
	return true;
}

vec Gain::get_row ( int index )
{
	return C.row(index);
}
		
vec Gain::get_column ( int index )
{
	return C.col(index);
}
		
void Gain::set_row ( int index,vec new_row_elements )
{
	//Its not implemented yet
}

void Gain::set_row ( int index,char* new_row_elements )
{
	vec new_row= vec(new_row_elements);
	set_row (index,new_row);
}

double Gain::at ( int index_x ,int index_y ){
	return C.at(index_x,index_y);
}