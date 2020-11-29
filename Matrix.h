#pragma once
#include "Vector.h"
using namespace std;

class Matrix
{
public:
   int N;               // Размер матрицы
   int M;               // Количество элементов в треугольнике
   
   vector<int> ig;      // Указатели начала строк
   vector<int> jg;      // Номера столбцов внедиагональных элементов

   vector<real> ggl;    // Верхний треугольник
   vector<real> ggu;    // Нижний треугольник

   vector<real> di;     // Диагональ

   Matrix(string path)
   {
      ifstream fin;
      fin.open(path + "kuslau.txt");
      fin >> N;
      fin.close();

      read_vector<int>(path + "ig.txt", ig, N + 1);
      M = ig[N];
      read_vector<int>(path + "jg.txt", jg, M);
      read_vector<real>(path + "ggl.txt", ggl, M);
      read_vector<real>(path + "ggu.txt", ggu, M);
      read_vector<real>(path + "di.txt", di, N);
   }

   Matrix(int _N, int _M)
   {
      N = _N;
      M = _M;

      ggl.resize(M);
      ggu.resize(M);
      jg.resize(M);
      di.resize(N);
      ig.resize(N + 1);
   }
   Matrix()
   {

   }

   // Получение диагональной факторизации матрицы
   Matrix diag_fact()
   {
      Matrix fact = Matrix(N, 0);

      fact.di = di;

      for (int i = 0; i < N + 1; i++)
         fact.ig[i] = 0;

      return fact;
   }

   // Получение неполного разложения Холецкого матрицы
   Matrix holec()
   {
      Matrix fact = Matrix(N, M);

      fact.ig = ig;
      fact.jg = jg;

      fact.di[0] = sqrt(di[0]);

      for (int i = 0; i < N; i++)
      {
         int prof_len = ig[i + 1] - ig[i];

         if (prof_len)
         {
            int i_in_gg = ig[i];
            if (!jg[i_in_gg])
               fact.ggl[i_in_gg] = fact.ggu[i_in_gg] = ggl[i_in_gg] / di[0];
         }
      }

      for (int i = 1; i < N; i++)
      {
         int prof_len = ig[i + 1] - ig[i];
         for (int k = 0; k < prof_len; k++)
         {
            int i_in_gg = ig[i] + k;
            int j = jg[i_in_gg];

            real sum = 0;

            for (int k_i = 0; k_i < k; k_i++)
            {
               int i_in_gg_2 = ig[i] + k_i;

               for (int k_j = 0; k_j < k; k_j++)
               {
                  int i_in_gg_3 = ig[j] + k_i;
                  if (jg[i_in_gg_3] == jg[i_in_gg_2])
                  {
                     sum += fact.ggl[i_in_gg_2] * fact.ggl[i_in_gg_3];
                     break;
                  }
               }
            }

            fact.ggl[i_in_gg] = fact.ggu[i_in_gg] = (ggl[i_in_gg] - sum) / fact.di[j];
         }

         real sum = 0;
         for (int k = 0; k < prof_len; k++)
         {
            int i_in_gg = ig[i] + k;
            sum += fact.ggl[i_in_gg] * fact.ggl[i_in_gg];
         }

         fact.di[i] = sqrt(abs(di[i] - sum));
      }
      return fact;
   }
};