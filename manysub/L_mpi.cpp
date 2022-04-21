#include <iostream>
#include <fstream>
#include <complex>
#include <boost/format.hpp>
#include <fmt/core.h>
#include <time.h>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/collectives.hpp>
#include "mfsw.hpp"
#include <algorithm>

using namespace std;
using namespace fmt;
namespace bmpi = boost::mpi;


struct param {
double J;
double K;
double h;
double hx;
double hy;
double hz;
double gamma;
double gamma_;
double T;
};

void xbond(int n, int i, int j, mfsw& ms, const param& p)
{
  int Nc = 3; // サイトあたりの平均場の数
  CPPL::dsymatrix Js(Nc);
  Js.zero();
  Js(0, 0) = p.J + 2. * p.K;  Js(0, 1) = p.gamma_;  Js(0, 2) = p.gamma_;
  Js(1, 0) = p.gamma_;       Js(1, 1) = p.J;       Js(1, 2) = p.gamma; 
  Js(2, 0) = p.gamma_;       Js(2, 1) = p.gamma;   Js(2, 2) = p.J;
  ms.set_J(n) = std::tuple<int, int, CPPL::dsymatrix>{i, j, Js}; // bond index 0 として、サイト0 と サイト1 の相互作用をセット
}
void ybond(int n, int i, int j, mfsw& ms, const param& p)
{
  int Nc = 3; // サイトあたりの平均場の数
  CPPL::dsymatrix Js(Nc);
  Js.zero();
  Js(0, 0) = p.J;       Js(0, 1) = p.gamma_;       Js(0, 2) = p.gamma;
  Js(1, 0) = p.gamma_;  Js(1, 1) = p.J + 2. * p.K;  Js(1, 2) = p.gamma_; 
  Js(2, 0) = p.gamma;   Js(2, 1) = p.gamma_;       Js(2, 2) = p.J;
  ms.set_J(n) = std::tuple<int, int, CPPL::dsymatrix>{i, j, Js}; // bond index 0 として、サイト0 と サイト1 の相互作用をセット
}
void zbond(int n, int i, int j, mfsw& ms, const param& p)
{
  int Nc = 3; // サイトあたりの平均場の数
  CPPL::dsymatrix Js(Nc);
  Js.zero();
  Js(0, 0) = p.J;       Js(0, 1) = p.gamma;   Js(0, 2) = p.gamma_;
  Js(1, 0) = p.gamma;   Js(1, 1) = p.J;       Js(1, 2) = p.gamma_; 
  Js(2, 0) = p.gamma_;  Js(2, 1) = p.gamma_;  Js(2, 2) = p.J + 2. * p.K;
  ms.set_J(n) = std::tuple<int, int, CPPL::dsymatrix>{i, j, Js}; // bond index 0 として、サイト0 と サイト1 の相互作用をセット
}

void site(int i, double rx, double ry, mfsw& ms)
{
  int Nc = 3;
  vector<double> spin(Nc);
  spin[0] = ms.mf_val( i * 3 );
  spin[1] = ms.mf_val( i * 3 + 1 );
  spin[2] = ms.mf_val( i * 3 + 2 );
  ms.set_site(i) = std::tuple< int, double, double, vector<double> >{i, rx, ry, spin}; //サイトの位置と平均場から得られたスピンの情報をセット
}


int main(int argc, char *argv[])
{
  bmpi::environment env(argc, argv);
  bmpi::communicator world;
  double pi = M_PI;
  double r3 = sqrt(3.);
//   struct timespec startTime, endTime;

//   if(argc!=5) exit(1);



  complex<double> im(0, 1);

  int Ns = 2; // 局所空間の次元
  CPPL::zhematrix Sxop(Ns);
  Sxop.zero();
  CPPL::zhematrix Syop(Ns);
  Syop.zero();
  CPPL::zhematrix Szop(Ns);
  Szop.zero();

  if (Ns == 2)
  {
    // cout << " S = 1/2" << endl;
    Sxop(0, 1) = 0.5;
    Syop(0, 1) = -0.5 * im;
    Szop(0, 0) = 0.5;
    Szop(1, 1) = -0.5;
  }
  else if (Ns == 4)
  {
    // cout << " S = 3/2" << endl;
    Sxop(0, 1) = r3 / 2.;
    Sxop(1, 2) = 1.;
    Sxop(2, 3) = r3 / 2.;
    Syop(0, 1) = -r3 * im / 2.;
    Syop(1, 2) = -1. * im;
    Syop(2, 3) = -r3 * im / 2.;
    Szop(0, 0) = 3./2.;
    Szop(1, 1) = 1./2.;
    Szop(2, 2) = -1./2.;
    Szop(3, 3) = -3./2.;
  }


  int Ni = atoi(argv[1]);
  vector<int> L_tmp(world.size());
  if (world.rank() == 0)
  {
    for (int  i = 0; i < L_tmp.size(); i++)
      std::cin >> L_tmp[i];
  }
  int L;
  bmpi::scatter(world, L_tmp, L, 0); 
  int M = 2 * L * L; //副格子 
  int N = M * 3;  // 平均場の数
  int Nc = 3; // サイトあたりの平均場の数

  int Nb = M * 3. / 2.; //ボンドの種類

  param p;
  double a = atof(argv[2]);  //K,Jの変数
//   p.h = atof(argv[3]); //印加磁場の絶対値
  p.T = 1e-8;            // 温度
  p.J = 0.;
  p.K = 1. / 2.;
  p.gamma = 0.5;
  p.gamma_ = -0.02;
  world.barrier();
  if (world.rank() == 0)
    cout << "J = " << p.J << " : K = " << 2 * p.K << " : Γ = " << p.gamma << " : Γ' = " << p.gamma_ << endl; 

//   std::string file = fmt::format("α ={}_gamma={}_gamma^={}.txt", a, p.gamma, p.gamma_);
//   std::string file = fmt::format("K=1_gamma={}_gamma^={}.txt", p.gamma, p.gamma_);
std::string file = fmt::format("tmp.txt");
  ofstream ofs(file);
  for (p.h = 2.; p.h > 1.; p.h -= 0.05)
  {
    mfsw ms(N, Nc, Nb, Ni, 1e-9);

    for (size_t i = 0; i < N; i += Nc)
    {
      ms.set_mat(i) = Sxop;
      ms.set_mat(i+1) = Syop;
      ms.set_mat(i+2) = Szop;
    }
  
    p.hx = p.h / sqrt(3.);
    p.hy = p.h / sqrt(3.);
    p.hz = p.h / sqrt(3.);

    int n_ = 0;
    for (int m = 1; m < 2 * L; m += 2)
    {
      for (int i = 0; i < L; i++)
      {
          xbond(n_, m*L+i, (m-1)*L+i, ms, p); n_++;
      }
      for (int i = 0; i < L; i++)
      {
        if (m*L+i == (m+1)*L-1)
        {
          ybond(n_, m*L+i, (m-1)*L, ms, p); n_++;
        }
        else
        {
          ybond(n_, m*L+i, (m-1)*L+(i+1), ms, p); n_++;
        }
      }
      for (int i = 0; i < L; i++)
      {
        if (m == 2*L-1)
        {
          zbond(n_, m*L+i, i, ms, p); n_++;
        }
        else
        {
          zbond(n_, m*L+i, (m+1)*L+i, ms, p); n_++;
        }
      }
    }
    
    for (size_t i = 0; i < N; i+=Nc)
    {
      ms.set_mat(i) = Sxop;
      ms.set_mat(i+1) = Syop;
      ms.set_mat(i+2) = Szop;
    }
    ms.set_T() = p.T;
    for(size_t i = 0; i < N; i += 3)
    {
      ms.set_H(i) = p.h / sqrt(3.);
      ms.set_H(i+1) = p.h / sqrt(3.);
      ms.set_H(i+2) = p.h / sqrt(3.);
    }
    std::string cand = fmt::format("../ct{}sub.txt", M);
    ms.exec_mf(cand);

    std::stringstream ss;
    ss << "###  I am thread " << world.rank() <<  ". I have size L = " << L << " : h = " << p.h << endl << ms.mf_out();
    for (int i = 0 ; i < world.size(); i++)
    {
      if (i == world.rank())
      {
        std::cout <<  ss.str() << std::endl;;
      }
    }
    world.barrier();
    // std::cout << ms.mf_out() << std::endl;
    double ene_tmp = ms.mf_ene();
    vector<double> ene(world.size());
    bmpi::gather(world, ene_tmp, ene, 0);
    if (world.rank() == 0)
    {
      double ene_min = *min_element(ene.begin(), ene.end());
      double L_out;
      for (int i = ene.size() - 1; i > 0 - 1; i--)
      {
        if ( abs( ene_min - ene[i] ) < 1e-5 )
          L_out = L_tmp[i];
      }
      ofs << p.h << " " << L_out << std::endl;
    }
  }


//   // ofstream ofs("fc.txt");
//   //ofstream ofs2("spin8sub.txt");
//   //ofstream ofs2(fn);

//   int d = 2;
//   vector <double> b1(d);
//   vector <double> b2(d);
//   vector <double> a1(d);
//   vector <double> a2(d);
//   vector <double> r1(d);
//   vector <double> r2(d);

//   a1[0] = 1. * (double)L;      a1[1] = 0. * (double)L;
//   a2[0] = 1./2. * (double)L;   a2[1] = r3/2. * (double)L;    //direct lattice vector

//   r1[0] = 1.;      r1[1] = 0.;
//   r2[0] = 1./2.;   r2[1] = r3/2.;    //direct lattice vector (2sub)

//   b1[0] = (4. * M_PI) / (r3 * (double)L) * r3/2.;  b1[1] = (4. * M_PI) / (r3 * (double)L) * -1./2.;
//   b2[0] = (4. * M_PI) / (r3 * (double)L) * 0;      b2[1] = (4. * M_PI) / (r3 * (double)L) * 1.;
//   int cntx = 0;
//   int cnty = 0;

//   int cnt = 0;
  // for (int n = 0; n < 2*L; n += 2)
  // {
  //   for (int i = 0; i < L; i++)
  //   {
  //     site(n * L + i, cos(pi/6.)/r3 + 1. * i + r2[0] * n/2, sin(pi/6.)/r3 + r2[1] * n/2, ms);
  //     cnt++;
  //   }
  // }
  // for (int n = 1; n < 2*L; n += 2)
  // {
  //   for (int i = 0; i < L; i++)
  //   {
  //     site(n * L + i, cos(11.*pi/6.)/r3 + 1. * i + r2[0] * (n-1)/2, sin(11.*pi/6.)/r3 + r2[1] * (n-1)/2, ms);
  //     cnt++;
  //   }
  // }
  // cout << "Total sites : " << cnt << endl;


  // for(int n = 0;n<4;n++)
  // {
  // clock_gettime(CLOCK_REALTIME, &startTime);
  //   for (double m1 = -3.; m1 <= 3.; m1 += 0.01)
  //   {
  //     cntx += 1;
  //     cnty = 0;
  //     for (double m2 = -3.; m2 <= 3.; m2 += 0.01)
  //     {

  //       cnty += 1;
  //       // double kx = m1 * b1[0] + m2 * b2[0];
  //       // double ky = m1 * b1[1] + m2 * b2[1];
  //       double kx = m1 * M_PI;
  //       double ky = m2 * M_PI;

  //       complex<double> gx = exp((0.5 * kx + 0.5 / sqrt(3.) * ky) * im);
  //       complex<double> gy = exp((-0.5 * kx + 0.5 / sqrt(3.) * ky) * im);
  //       complex<double> gz = exp((-1. / sqrt(3.) * ky) * im);
  //       complex<double> dkx_gx = 0.5 * im * gx;
  //       complex<double> dkx_gy = -0.5 * im * gy;
  //       complex<double> dkx_gz = 0.;
  //       complex<double> dky_gx = 0.5/sqrt(3) * im * gx;
  //       complex<double> dky_gy = 0.5/sqrt(3) * im * gy;
  //       complex<double> dky_gz = -1./sqrt(3) * im * gz;
  //       vector<complex<double>> gamma = {gx,gx, gy,gy, gz,gz};
  //       vector<complex<double>> dkx_gamma = {dkx_gx,dkx_gx, dkx_gy,dkx_gy, dkx_gz,dkx_gz};
  //       vector<complex<double>> dky_gamma = {dky_gx,dky_gx, dky_gy,dky_gy, dky_gz,dky_gz};
  //       ofs << kx << " " << ky << " " << ms.fc(a1, a2, kx, ky) << endl;
  //     }
  //   }
  // clock_gettime(CLOCK_REALTIME, &endTime);

  // printf("elapsed time = ");
  //     if (endTime.tv_nsec < startTime.tv_nsec) {
  //     printf("%5ld.%09ld", endTime.tv_sec - startTime.tv_sec - 1,
  //             endTime.tv_nsec + (long int)1.0e+9 - startTime.tv_nsec);
  //     } else {
  //     printf("%5ld.%09ld", endTime.tv_sec - startTime.tv_sec,
  //             endTime.tv_nsec - startTime.tv_nsec);
  //     }
  //     printf("(sec)\n");

  // }
}
