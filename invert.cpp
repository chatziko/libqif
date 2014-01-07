#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char** argv) {
  mat A;
  A.load(argv[1], arma_ascii);
  mat i = A.i();
  
  cout << A.n_cols;
  i.save(argv[2], arma_ascii);
  
  return 0;
}
