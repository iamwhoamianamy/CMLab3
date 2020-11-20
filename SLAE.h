#pragma once
#include <fstream>
#include <vector>
using namespace std;

typedef double real;

class SLAE
{
public:
   int N;            // Размер матрицы
   int M;            // Количество ненулевых элементов в треугольнике
   int maxiter;      // Максимальное количество итераций
   real eps;         // Велечина требуемой относительной невязки

   vector<int> ig;   // Указатели начала строк
   vector<int> jg;   // Номера столбцов внедиагональных элементов
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

      read_vector<int>(path + "ig.txt", ig, N);
      M = ig[N];
      read_vector<int>(path + "jg.txt", jg, M);
      read_vector<real>(path + "ggl.txt", ggl, M);
      read_vector<real>(path + "ggu.txt", ggu, M);
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
      res.resize(N);

      for (int i = 0; i < N; i++)
      {

      }
   }
};

