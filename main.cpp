#include <iostream>
#include "SLAE.h"

using namespace std;

int main()
{
   SLAE slae = SLAE("test3/");

   vector<real> x0(slae.N, 0);
   vector<real> res(slae.N, 0);

   slae.conj_grad_method(x0, res);

   for (size_t i = 0; i < slae.N; i++)
   {
      cout << res[i] << " ";
   }
}