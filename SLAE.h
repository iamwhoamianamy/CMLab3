#pragma once
#include "Vector.h"
using namespace std;

typedef double real;

class SLAE
{
public:
   int N;               // ������ �������
   int M;               // ���������� ��������� � ������������
                        
   int maxiter;         // ������������ ���������� ��������
   real eps;            // �������� ��������� ������������� �������
                        
   vector<int> ig;      // ��������� ������ �����
   vector<int> jg;      // ������ �������� ��������������� ���������
                       
   vector<real> ggl;    // ������� �����������
   vector<real> ggu;    // ������ �����������
                       
   vector<real> di;     // ���������
   vector<real> pr;     // ������ ������ �����
                       
   vector<real> t;      // ��������������� ������ ��� ���
   vector<real> rk;     // ������ ������� �� ������� �������� ���
   vector<real> rk1;    // ������ ������� �� ����. �������� ���
   vector<real> zk1;    // ������ ������ �� ����. �������� ���
   vector<real> AtAzk1; // ��������������� ������ ��� ���

   // ����������� ������ SLAE
   SLAE(string path)
   {
      ifstream fin;

      fin.open(path + "kuslau.txt");
      fin >> N >> maxiter >> eps;
      fin.close();

      read_vector<int>(path + "ig.txt", ig, N + 1);
      M = ig[N];
      read_vector<int>(path + "jg.txt", jg, M);
      read_vector<real>(path + "ggl.txt", ggl, M);
      read_vector<real>(path + "ggu.txt", ggu, M);
      read_vector<real>(path + "di.txt", di, N);
      read_vector<real>(path + "pr.txt", pr, N);

      t.resize(N);
      rk.resize(N);
      rk1.resize(N);
      zk1.resize(N);
      AtAzk1.resize(N);
   }

   // ������� ��������� ������� �� ������ vec, ��������� � res
   void matrix_vector_mult(const vector<real>& vec, vector<real>& res,
                           const vector<real>& bot_tr, const vector<real>& top_tr)
   {
      for (int i = 0; i < N; i++)
         res[i] = 0;

      for (int i = 0; i < N; i++)
      {
         res[i] += vec[i] * di[i];

         int prof_len = ig[i + 1] - ig[i];
         for (int k = 0; k < prof_len; k++)
         {
            int index_in_prof = ig[i] + k;
            int j = jg[index_in_prof];
            res[i] += vec[j] * bot_tr[index_in_prof];
            res[j] += vec[i] * top_tr[index_in_prof];
         }
      }
   }

   // ����� ����������� ����������, ���������� ���������� ��������
   int conj_grad_method(vector<real> xk1, vector<real>& res)
   {
      for (int i = 0; i < N; i++)
         res[i] = 0;

      matrix_vector_mult(xk1, t, ggl, ggu);  // t = A * x0
      matrix_vector_mult(pr - t, rk1, ggu, ggl);   // r0 = AT(f - A * x0)
      zk1 = rk1;

      int k = 1;
      while (k < maxiter)
      {
         matrix_vector_mult(zk1, t, ggl, ggu);     // t = A * zk-1

         matrix_vector_mult(t, AtAzk1, ggu, ggl); // AtAzk1 = At * A * zk-1
         real ak = (rk1 * rk1) / (AtAzk1 * zk1);
         xk1 = xk1 + ak * zk1;
         rk = rk1 - ak * AtAzk1;
         real bk = (rk * rk) / (rk1 * rk1);
         zk1 = rk + bk * zk1;
         rk1 = rk;
         real disc = norm(rk) / norm(pr);

         if (disc < eps)
            break;
         else
            k++;
      }

      res = xk1;
      return k;
   }
};

