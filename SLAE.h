#pragma once
#include <fstream>
#include <vector>
using namespace std;

typedef double real;

class SLAE
{
public:
   int N;            // ������ �������
   int M;            // ���������� ��������� ��������� � ������������
   int maxiter;      // ������������ ���������� ��������
   real eps;         // �������� ��������� ������������� �������

   vector<int> ig;   // ��������� ������ �����
   vector<int> jg;   // ������ �������� ��������������� ���������
   vector<real> ggl; // ������� �����������
   vector<real> ggu; // ������ �����������
   vector<real> di;  // ���������
   vector<real> pr;  // ������ ������ �����
   
   // ����������� ������ SLAE
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

   // ���� ������� vec ����������� n �� ����� � ������ file_name
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

   // ������� ��������� ������� �� ������ vec, ��������� � res
   void matrix_vector_mult(const vector<real>& vec, vector<real>& res)
   {
      res.resize(N);

      for (int i = 0; i < N; i++)
      {

      }
   }
};

