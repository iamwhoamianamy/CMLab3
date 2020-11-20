#include <iostream>
#include "SLAE.h"

using namespace std;

int main()
{
   SLAE slae = SLAE("test1/");

   vector<real> res(4);

   slae.matrix_vector_mult(slae.pr, res);
}