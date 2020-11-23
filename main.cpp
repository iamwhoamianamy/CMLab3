#include <iostream>
#include "SLAE.h"

using namespace std;

int main()
{
   SLAE slae = SLAE("test1/");

   vector<real> res(4, 0);

   slae.conj_grad_method(res);

   for (size_t i = 0; i < 4; i++)
   {
      cout << res[i] << " ";
   }
}