#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <boost/format.hpp>

#include "mfsw.hpp"

using namespace std;


struct param {
double J;
double K;
double hx;
double hy;
double hz;
};

void xbond(int n, int i, int j, mfsw& ms, const param& p)
{
  int Nc = 3; // サイトあたりの平均場の数
  CPPL::dsymatrix Js(Nc);
  Js.zero();
  Js(0, 0) = p.J + 2 * p.K; // 相互作用行列をセット
  Js(1, 1) = p.J;
  Js(2, 2) = p.J;
  ms.set_J(n) = std::tuple<int, int, CPPL::dsymatrix>{i, j, Js}; // bond index 0 として、サイト0 と サイト1 の相互作用をセット
}
void ybond(int n, int i, int j, mfsw& ms, const param& p)
{
  int Nc = 3; // サイトあたりの平均場の数
  CPPL::dsymatrix Js(Nc);
  Js.zero();
  Js(0, 0) = p.J; // 相互作用行列をセット
  Js(1, 1) = p.J + 2 * p.K;
  Js(2, 2) = p.J;
  ms.set_J(n) = std::tuple<int, int, CPPL::dsymatrix>{i, j, Js}; // bond index 0 として、サイト0 と サイト1 の相互作用をセット
}
void zbond(int n, int i, int j, mfsw& ms, const param& p)
{
  int Nc = 3; // サイトあたりの平均場の数
  CPPL::dsymatrix Js(Nc);
  Js.zero();
  Js(0, 0) = p.J; // 相互作用行列をセット
  Js(1, 1) = p.J;
  Js(2, 2) = p.J + 2 * p.K;
  ms.set_J(n) = std::tuple<int, int, CPPL::dsymatrix>{i, j, Js}; // bond index 0 として、サイト0 と サイト1 の相互作用をセット
}



int main(int argc, char *argv[])
{


  complex<double> im(0, 1);

  int Ns = 2; // 局所空間の次元
  CPPL::zhematrix Sxop(Ns);
  Sxop.zero();
  CPPL::zhematrix Syop(Ns);
  Syop.zero();
  CPPL::zhematrix Szop(Ns);
  Szop.zero();

  Sxop(0, 1) = 0.5;
  Syop(0, 1) = -0.5 * im;
  Szop(0, 0) = 0.5;
  Szop(1, 1) = -0.5;

  int N = 24*3;  // 平均場の数
  int SL = 24; //sublatticeの数
  int Nc = 3; // サイトあたりの平均場の数

  int Nb = 24*3/2; //ボンドの種類

  mfsw ms(N, Nc, Nb, 1000, 1e-9);

  for (size_t i = 0; i < N; i+=Nc)
  {
  ms.set_mat(i) = Sxop;
  ms.set_mat(i+1) = Syop;
  ms.set_mat(i+2) = Szop;
  }

  param p;


  double a = atof(argv[1]);  //K,Jの変数
  double theta = atof(argv[2]); //磁場の角度依存
  double phi = atof(argv[3]); //磁場の角度依存 phi
  double h = atof(argv[4]); //印加磁場の絶対値
  p.hx = atof(argv[5]);
  p.hy = atof(argv[6]);
  p.hz = atof(argv[7]);
  
  p.J = cos(a);
  p.K = sin(a);

  int n=0;
  xbond(n, 0, 1, ms, p); n++;
  xbond(n, 2, 22, ms, p); n++;
  xbond(n, 4, 3, ms, p); n++;
  xbond(n, 21, 5, ms, p); n++;
  xbond(n, 6, 12, ms, p); n++;
  xbond(n, 7, 8, ms, p); n++;
  xbond(n, 9, 23, ms, p); n++;
  xbond(n, 11, 10, ms, p); n++;
  xbond(n, 13, 14, ms, p); n++;
  xbond(n, 15, 19, ms, p); n++;
  xbond(n, 17, 16, ms, p); n++;
  xbond(n, 20, 18, ms, p); n++;
  ybond(n, 0, 5, ms, p); n++;
  ybond(n, 20, 1, ms, p); n++;
  ybond(n, 2, 3, ms, p); n++;
  ybond(n, 4, 23, ms, p); n++;
  ybond(n, 6, 14, ms, p); n++;
  ybond(n, 7, 12, ms, p); n++;
  ybond(n, 21, 8, ms, p); n++;
  ybond(n, 9, 10, ms, p); n++;
  ybond(n, 11, 19, ms, p); n++;
  ybond(n, 13, 18, ms, p); n++;
  ybond(n, 15, 16, ms, p); n++;
  ybond(n, 17, 22, ms, p); n++;
  zbond(n, 0, 19, ms, p); n++;
  zbond(n, 2, 1, ms, p); n++;
  zbond(n, 6, 3, ms, p); n++;
  zbond(n, 4, 5, ms, p); n++;
  zbond(n, 7, 22, ms, p); n++;
  zbond(n, 9, 8, ms, p); n++;
  zbond(n, 20, 10, ms, p); n++;
  zbond(n, 11, 12, ms, p); n++;
  zbond(n, 13, 23, ms, p); n++;
  zbond(n, 15, 14, ms, p); n++;
  zbond(n, 21, 16, ms, p); n++;
  zbond(n, 17, 18, ms, p); n++;
  
  

  
  for(size_t i = 0; i < N; i += 3)
  {
    //ms.set_H(i) = p.h/sqrt(3.);
    ms.set_H(i) = p.hx;
    ms.set_H(i+1) = p.hy;
    ms.set_H(i+2) = p.hz;
  }


  //ofstream ofs1("mfhd24sub.txt",ios::app);
  // 初期値のinputファイル(行数分だけ回る。列数はNと合致させる)を用いて平均場近似を実行
  ms.exec_mf("cand24sub_hd.txt");  //forced ferr



  cout << "# ";
  cout << ms.mf_out() << endl;

}
