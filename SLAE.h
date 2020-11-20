#pragma once
#include <fstream>
#include <vector>
using namespace std;

typedef double real;

class SLAE
{
public:
   int N;            // Размер матрицы

   int Ml;           // Количество ненулевых элементов в нижнем треугольнике
   int Mu;           // Количество ненулевых элементов в верхнего треугольнике
  
   int maxiter;      // Максимальное количество итераций
   real eps;         // Велечина требуемой относительной невязки

   vector<int> igl;  // Указатели начала строк нижнего треугольника
   vector<int> igu;  // Указатели начала строк верхнего треугольника

   vector<int> jgl;  // Номера столбцов внедиагональных элементов нижнего треуголника
   vector<int> jgu;  // Номера столбцов внедиагональных элементов верхнего треуголника
   
   vector<real> ggl; // Верхний треугольник
   vector<real> ggu; // Нижний треугольник
   
   vector<real> di;  // Диагональ
   vector<real> pr;  // Вектор правой части
   
   // Конструктор класса SLAE
   SLAE(string path)
   {
      ifstream fin;

      fin.open(path + "kuslau.txt");
      fin >> N >> maxiter >> eps;
      fin.close();

      read_vector<int>(path + "igl.txt", igl, N + 1);
      Ml = igl[N];
      read_vector<int>(path + "igu.txt", igu, N + 1);
      Mu = igu[N];

      read_vector<int>(path + "jgl.txt", jgl, Ml);
      read_vector<int>(path + "jgu.txt", jgu, Mu);

      read_vector<real>(path + "ggl.txt", ggl, Ml);
      read_vector<real>(path + "ggu.txt", ggu, Mu);

      read_vector<real>(path + "di.txt", di, N);
      read_vector<real>(path + "pr.txt", pr, N);
   }

   // Ввод вектора vec размерности n из файла с именем file_name
   template<class type>
   void read_vector(string file_name, vector<type>& vec, int n)
   {
      vec.resize(n);
      ifstream fin;
      fin.open(file_name);
      for (int i = 0; i < n; i++)
         fin >> vec[i];

      fin.close();
   }

   // Функция умножения матрицы на вектор vec, результат в res
   void matrix_vector_mult(const vector<real>& vec, vector<real>& res)
   {
      for (int i = 0; i < N; i++)
      {
         // Диагональ
         res[i] += vec[i] * di[i];

         // Нижний треугольник
         int prof_len = igl[i + 1] - igl[i];
         for (int k = 0; k < prof_len; k++)
         {
            int index_in_prof = igl[i] + k;
            int j = jgl[index_in_prof];
            res[i] += vec[j] * ggl[index_in_prof];
         }

         // Верхний треугольник
         prof_len = igu[i + 1] - igu[i];
         for (int k = 0; k < prof_len; k++)
         {
            int i_in_prof = igu[i] + k;
            int j = jgu[i_in_prof];
            res[i] += vec[j] * ggu[i_in_prof];
         }
      }
   }
};

