#pragma once
#include "Matrix.h"
using namespace std;

class SLAE
{
public:                    
   int maxiter;         // Максимальное количество итераций
   real eps;            // Велечина требуемой относительной невязки

   Matrix mat;          // Матрица

   vector<real> pr;     // Вектор правой части

   vector<real> t;      // Вспомогательный вектор для МСГ
   vector<real> tt;     // Вспомогательный вектор для МСГ
   vector<real> rk1;    // Вектор невязки на перд. итерации МСГ
   vector<real> zk1;    // Вектор спуска на пред. итерации МСГ
   vector<real> AtAzk1; // Вспомогательный вектор для МСГ

   SLAE()
   {

   }

   SLAE(string path)
   {
      ifstream fin;

      int N;

      fin.open(path + "kuslau.txt");
      fin >> N >> maxiter >> eps;
      fin.close();

      mat = Matrix(path);

      pr.resize(N);

      t.resize(N);
      tt.resize(N);
      rk1.resize(N);
      zk1.resize(N);
      AtAzk1.resize(N);
   }

   SLAE(int N, int _maxiter, real _eps, const Matrix& _mat)
   {
      maxiter = _maxiter;
      eps = _eps;

      mat = Matrix(_mat);

      pr.resize(N);

      t.resize(N);
      rk1.resize(N);
      zk1.resize(N);
      AtAzk1.resize(N);
   }

   // Функция умножения матрицы на вектор vec, результат в res
   void matrix_vector_mult(const vector<real>& vec, vector<real>& res,
                           const vector<real>& bot_tr, const vector<real>& top_tr)
   {
      for (int i = 0; i < mat.N; i++)
         res[i] = 0;

      for (int i = 0; i < mat.N; i++)
      {
         res[i] += vec[i] * mat.di[i];

         int prof_len = mat.ig[i + 1] - mat.ig[i];
         for (int k = 0; k < prof_len; k++)
         {
            int i_in_gg = mat.ig[i] + k;
            int j = mat.jg[i_in_gg];
            res[i] += vec[j] * bot_tr[i_in_gg];
            res[j] += vec[i] * top_tr[i_in_gg];
         }
      }
   }


   // Метод сопряженных градиентов, возвращает количество итераций
   int conj_grad_method(vector<real>& xk1, vector<real>& res)
   {
      for (int i = 0; i < mat.N; i++)
         res[i] = xk1[i] = 0;

      matrix_vector_mult(xk1, t, mat.ggl, mat.ggu);       // t = A * x0
      matrix_vector_mult(pr - t, rk1, mat.ggu, mat.ggl);  // r0 = AT(f - A * x0)
      zk1 = rk1;

      int k = 1;
      while (k < maxiter)
      {
         matrix_vector_mult(zk1, t, mat.ggl, mat.ggu);    // t = A * zk-1

         matrix_vector_mult(t, AtAzk1, mat.ggu, mat.ggl); // AtAzk1 = At * A * zk-1
         real ak = (rk1 * rk1) / (AtAzk1 * zk1);
         xk1 = xk1 + ak * zk1;
         real bk = rk1 * rk1;
         rk1 = rk1 - ak * AtAzk1;
         bk = (rk1 * rk1) / bk;
         zk1 = rk1 + bk * zk1;

         real disc = norm(rk1) / norm(pr);

         if (disc < eps)
            break;
         else
            k++;
      }

      res = xk1;
      return k;
   }

   // Метод сопряженных градиентов с предобусловденной неполной
   // факторизацией матрицей, возвращает количество итераций
   int conj_grad_pred_method(vector<real>& xk1, vector<real>& res, SLAE& fac_slae)
   {
      for (int i = 0; i < mat.N; i++)
         res[i] = xk1[i] = 0;

      matrix_vector_mult(xk1, t, mat.ggl, mat.ggu);       // t = A * x0
      matrix_vector_mult(pr - t, rk1, mat.ggu, mat.ggl);  // r0 = AT(f - A * x0)

      // Решаем z0 = M-1 * r0
      fac_slae.pr = rk1;
      fac_slae.conj_grad_method(t, zk1);
      
      int k = 1;
      while (k < maxiter)
      {
         // Решаем tt = M-1 * rk-1
         fac_slae.pr = rk1;
         fac_slae.conj_grad_method(t, tt);

         matrix_vector_mult(zk1, t, mat.ggl, mat.ggu);    // t = A * zk-1
         matrix_vector_mult(t, AtAzk1, mat.ggu, mat.ggl); // AtAzk1 = At * A * zk-1

         real ak = (tt * rk1) / (AtAzk1 * zk1);
         xk1 = xk1 + ak * zk1;
         real bk = tt * rk1;
         rk1 = rk1 - ak * AtAzk1;

         // Решаем tt = M-1 * rk
         fac_slae.pr = rk1;
         fac_slae.conj_grad_method(t, tt);

         bk = (tt * rk1) / bk;
         zk1 = tt + bk * zk1;

         real disc = norm(rk1) / norm(pr);

         if (disc < eps)
            break;
         else
            k++;
      }

      res = xk1;
      return k;
   }
};

