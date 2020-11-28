#include <iostream>
#include <string>
#include <direct.h>
#include "SLAE.h"

using namespace std;

// Генерация матрицы Гильберта размерности n
void gen_hilb(int n)
{
   string path0 = "tests/hilb/";

   ofstream fout;
   fout.open(path0 + "kuslau.txt");

   fout << n << " " << "10000" << " " << "10e-14";

   fout.close();

   int ggl_len = n * (n - 1) / 2;

   vector<real> di(n);
   vector<int> ig(n + 1), jg(ggl_len);
   vector<real> ggl(ggl_len);

   ig[0] = ig[1] = 0;

   di[0] = 1;

   for (int i = 1; i < n; i++)
   {
      int i0 = ig[i + 0];
      int i1 = ig[i + 1] = ig[i + 0] + i;

      for (int j = 0, k = i0; j < i; j++, k++)
      {
         ggl[k] = 1.0 / ((real(i) + real(j) + 1.0));
         jg[k] = j;
      }

      di[i] = 1 / (real(i) * 2 + 1);
   }

   print_vector(path0 + "ggl.txt", ggl);
   print_vector(path0 + "ggu.txt", ggl);

   print_vector(path0 + "di.txt", di);
   print_vector(path0 + "ig.txt", ig);
   print_vector(path0 + "jg.txt", jg);
}

// Формирование отчета путем решение СЛАУ всеми методами
void report(string test, bool calc_pr)
{
   // Считываение матрицы из соответствующей папки
   string path0 = "tests/" + test + "/";  // Путь к папке с матрицей
   SLAE slae = SLAE(path0);               // Ввод матрицы из файлов
   int n = slae.mat.N;

   vector<real> x(n);
   vector<real> res(n);

   // Создание вектора x = 1,2,3,...,n , умножение матрицы на него,
   // записсcь результата в файл "pr.txt"
   if (calc_pr)
   {
      for (int i = 0; i < n; i++)
         x[i] = i + 1;

      slae.matrix_vector_mult(x, res, slae.mat.ggl, slae.mat.ggu);
      print_vector(path0 + "pr.txt", res);
   }
   
   // Чтение вектора правой части из файла
   read_vector(path0 + "pr.txt", slae.pr, n);

   // Формирование отчета в файл "report_.txt" в папке reports
   ofstream fout;
   fout.open("reports/report_" + test + ".txt");

   int k_iter;

   for (int i = 0; i < 3; i++)
   {
      // Выбор метода решения
      switch (i)
      {
      case 0:
      {
         // Решение методом сопряженных градиентов
         k_iter = slae.conj_grad_method(x, res);

         fout << "МСГ" << endl;
         break;
      }
      case 1:
      {
         // Решение методом сопряженных градиентов c неполной
         // диагональной факторизацией
         Matrix fac_mat = slae.mat.diag_fact();
         SLAE fac_slae = SLAE(slae.mat.N, slae.maxiter, slae.eps, fac_mat);
         k_iter = slae.conj_grad_pred_method(x, res, fac_slae);

         fout << "МСГ c неполной диагональной факторизацией" << endl;
         break;
      }

      case 2:
      {
         // Решение методом сопряженных градиентов c неполной
         // факторизацией Холецкого
         Matrix fac_mat = slae.mat.holec();
         SLAE fac_slae = SLAE(slae.mat.N, slae.maxiter, slae.eps, fac_mat);
         k_iter = slae.conj_grad_pred_method(x, res, fac_slae);

         fout << "МСГ c неполной факторизацией Холецкого" << endl;
         break;
      }
      }

      // Вывод информации в файл "report.txt"
      fout << k_iter << endl; // Количество итераций
      fout << setprecision(14);

      // Вывод вектора х
      for (int j = 0; j < n; j++)
         fout << res[j] << endl;

      if (calc_pr)
      {
         fout << endl;

         // Вывод погрешности
         for (int j = 0; j < n; j++)
            fout << real(j) + 1.0 - res[j] << endl;
      }
      fout << endl;
   }

   fout.close();
}

int main()
{
   // Работа тестовыми матрицами маленькой размерности
   report("small", 1);
   report("small_sym", 1);

   // Работа с матрицами Гильберта разной размерности
   for (int i = 0; i < 1; i++)
   {
      int n;
      cout << "Enter the size of " << i + 1 << " Hilbert matrix: ";
      cin >> n;
      gen_hilb(n);
      report("hilb", 1);
   }

   // Работа с матрицами с диагональным преобладанием
   report("diagdom", 1);
   report("diagdommin", 1);

   // Работа с большими матрицами
   report("big1", 0);
   //report("big2", 0);
}