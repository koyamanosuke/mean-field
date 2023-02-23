#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include "cpplapack/cpplapack.h"

using namespace std;


#include "mfsw.hpp"
#include "makedata.hpp"

using exch_type = CPPL::dgematrix;
struct param 
{
  vector<tuple<vector<double>, exch_type>> J;
  vector<vector<double>> H;
  double T;
};

void bond(int n, int Nc, int i, int j, const vector<double> &xs, mfsw<exch_type>& ms, 
          param &p)
{
  int m;
  exch_type Js;
  for (int i = 0; i < p.J.size(); i++)
  {
    vector<double> x = std::get<0>(p.J[i]);
    double sum = 0.;
    for (int k = 0; k < xs.size(); k++)
      sum += abs(x[k] - xs[k]);
    double eps = 1e-4;
    if (sum < eps)
    {
      m = i;
      Js = std::get<1>(p.J[m]);
      break;
    }
    else 
    {
      Js.resize(Nc, Nc); Js.zero();
    }
  }
  ms.set_J(n) = std::tuple<int, int, exch_type>{i, j, Js};
}


int main(int argc, char *argv[])
{
  if(argc!=1) 
    exit(1);

  std::string fn_lattice = "data/input/spin_wave_230203/cubic2subbcc_lattice.dat";
  std::string fn_iij = "data/input/spin_wave_230203/cubic2subbcc_iij_1_1.dat";
  std::string fn_base = "data/input/spin_wave_230203/cubic2subbcc_base_1_1.dat";

  makedata<exch_type> md(fn_lattice, fn_iij, fn_base);

  complex<double> im(0, 1);

  int Ns = md.make_Ns();       // 局所状態の数  (generator の次元)
  int Nc = md.make_Nc();       // サイトあたりの平均場の数
  int M = 2;                   // 副格子の数
  int N = M * Nc;              // 平均場の数
  int Z = 8;                   // １つの注目サイトと相互作用するサイトの合計（最近接間でのみ相互作用する bcc なら　8）
  int Nb = M * Z/2;            // ユニットセルに含まれるボンドの総数
  int Nr = 5;                  // 初期値を（ランダムに）入力する回数
  int Ni = 100000;             // ループ回数



  ////*generator*////
  vector<CPPL::zhematrix> Oop(Nc);   // generator
  for (int i = 0; i < Oop.size(); i++)
    Oop[i] = md.make_generator()[i];


  param p;
  md.make_J(p.J);
  md.make_H(p.H, M);

  mfsw<exch_type> ms(N, Nc, Nb, Ni, Nr, 1e-8);

  for (size_t i = 0; i < N; i += Nc)
  {
    for (size_t j = 0; j < Nc; j++)
      ms.set_mat(i + j) = Oop[j];
  }

  for(size_t i = 0; i < N; i += Nc)
  {
    for (size_t j = 0; j < Nc; j++)
    {
      ms.set_H(i + j) = p.H[i / Nc][j];
    }
  }


  vector<vector<double>> xs(Nb);   // ユニットセル内の独立なボンド間の相対座標

  /* bond input */
  {
    ms.set_bond(xs, Nb, md);       // 元データと同じ副格子構造を仮定していれば md.set_bond を呼び出すだけで良い
  }

  ms.exec_mf();                    // Nr 回のランダムな初期値で解を探す
  cout << "# ";
  cout << ms.mf_out() << endl;     // 平均場解の出力（Nr回実行した内の最安定解のみ出力される）




  p.T = 1e-4;                      // System の温度
  ms.set_T() = p.T;
  std::function<complex<double>(vector<double>, vector<double>)> g_;
  g_ = [&im](vector<double> x, vector<double> k)
  { 
    complex<double> pd(1);
    for (int i = 0; i < x.size(); i++)
      pd = pd * exp(im * x[i] * k[i]);

    return pd;
  };

  vector<std::function<complex<double>(vector<double>, vector<double>)>> g(Nb);
  for (size_t i = 0; i < Nb; i++)
    g[i] = g_;

  vector<vector<double>> a(3);         // 基本並進ベクトル
  md.make_unitvec(a);

  std::string fnsw = "data/output/spec.txt";
  std::ofstream sw(fnsw);

  {    
    double x = 0.;
    double t = 0.;

    /* (-2π,-2π,-2π) --- (2π,2π,2π) line */
    for (x = -2.; x < 2.001; x += 0.005)
    {
      vector<double> k(3);
      k[0] = x;
      k[1] = x;
      k[2] = x;
      k[0] *= M_PI;
      k[1] *= M_PI;
      k[2] *= M_PI;
      vector<complex<double>> gamma(g.size());
      for (int i = 0; i < Nb; i++)
        gamma[i] = g[i](xs[i], k);
      sw << x << " " << ms.exec_sw_out(gamma);   // スピン波エネルギーおよび動的構造因子の書き込み
      sw << endl;
    }
    t = x;
  }

  std::string fnmg =  "data/output/swmag.txt";
  std::ofstream mg(fnmg);
  {
    double dk = 0.05;
    mg << " " << ms.exec_sw_mag_out(xs, g, a, dk) << endl;  // SW モーメントの書き込み
  }

}
