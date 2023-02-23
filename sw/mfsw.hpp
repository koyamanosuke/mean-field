#ifndef MFSW_HPP
#define MFSW_HPP

//////////////////////////////

// v 1.0
// 2023 / 2 / 1
//////////////////////////////

#include <iostream>
#include <complex>
#include <boost/format.hpp>
#include <iterator>
#include <random>

#include "bogoliubov.hpp"
#include "cpplapack/cpplapack.h"
#include "makedata.hpp"

using namespace std;

class RND : public std::mt19937{
private:
  unsigned int v_max;
public:
  RND(unsigned int n) : std::mt19937(n),v_max(std::mt19937::max()){};
  RND() : v_max(std::mt19937::max()){};
  int N(int n){ return (int)((unsigned)operator()()%n); }
  double d(){ return (double)operator()()/((double)v_max); }
  double d(double max){ return (double)operator()()/((double)v_max) * max; }
};

class mf
{
private:
  double val;
  double MIX;
  CPPL::zhematrix mat;

public:
  CPPL::zhematrix &assign(double MIX_)
  {
    MIX = MIX_;
    return mat;
  }
  void set_init(const double init_)
  {
    val = init_;
  }
  CPPL::zhematrix op()
  {
    return mat;
  }

  // 平均場の元での期待値
  double calc(CPPL::zhematrix &mat, CPPL::zcovector &v)
  {
    CPPL::zcovector vt = mat * v;
    complex<double> sum(0., 0.);
    for (int i = 0; i < mat.m; i++)
        sum += conj(v(i)) * vt(i);
    return real(sum);
  }

  double calc(const CPPL::zcovector &v)
  {
    CPPL::zcovector vt = mat * v;
    complex<double> sum(0, 0);
    for (int i = 0; i < v.l; i++)
      sum += conj(v(i)) * vt(i);
    double val_old = val;
    val = MIX * std::real(sum) + (1 - MIX) * val_old;
    return fabs(val - val_old);
  }
  double ave() { return val; };
};

template<typename S>
class mfsw
{

private:
  const int Nitr;
  const double MIX;
  const double EPS;
  const int Nr;
  using exch_type = S;

  int N;      // 平均場の数
  int M;      // 副格子の数
  int Nc;     // サイトあたりの平均場の数
  int Nb;     // ボンドの種類
  double T;   // 温度

  vector<tuple<int, int, exch_type>> bond; // i,jサイトを繋ぐ相互作用行列
  vector<double> hs;                             // iサイトに働く外場


  vector<mf> MF;

  vector<double> ave;
  vector<vector<double>> w;                 // 局所固有エネルギー
  vector<vector<CPPL::zcovector>> v;        // 局所状態

  vector<double> ene;
  double ene_min;
  


  CPPL::zgematrix sw_U;
  CPPL::zgematrix sw_V;
  CPPL::zgematrix sw_Jinv;
  vector<complex<double>> Drate;
  vector<complex<double>> Dself;
  vector<vector<vector<complex<double > > > > Vbint;
  vector<vector<vector<complex<double > > > > Wbint;
  vector<vector<vector<complex<double > > > > Cbint;

  template <typename T>
  void read_file(std::string filename, std::vector<std::vector<T>> &x)
  {
    x.clear();
    std::string line;
    std::ifstream ifs(filename.c_str());
    if (!ifs)
    {
      std::cerr << "Can't open the file : " << filename << std::endl;
      std::exit(1);
    }
    int i = 0;
    while (!ifs.eof())
    {
      std::vector<T> temp;
      getline(ifs, line);
      std::istringstream is(line);
      std::copy(std::istream_iterator<T>(is), std::istream_iterator<T>(), std::back_inserter(temp));
      if (temp.size())
      {
        x.resize(i + 1);
        x[i] = temp;
        i++;
      }
    }
  }

  double dot(const vector<double> &x, const vector<double> &y)
  {
    double sum = 0.;
    for (int i = 0; i < x.size(); i++)
      sum += x[i] * y[i];
    return sum;
  }
  vector<double> cross(const vector<double> &x, const vector<double> &y)
  {
    int N = x.size();
    std::vector<double> c(N);
    c[0] = x[1] * y[2] - x[2] * y[1];
    c[1] = x[2] * y[0] - x[0] * y[2];
    c[2] = x[0] * y[1] - x[1] * y[0];
    return c;
  }
  void reciprocal(vector<vector<double>> &a, vector<vector<double>> &b)
  {
    int Nd = a.size();
    if (Nd == 1)
    {
      a.resize(3);
      a[1] = {0, 1, 0};
      a[2] = {0, 0, 1};
    }
    if (Nd == 2)
    {
      a.resize(3);
      a[2] = {0, 0, 1};
    }

    b.resize(3);
    b[0] = cross(a[1], a[2]);
    b[1] = cross(a[2], a[0]);
    b[2] = cross(a[0], a[1]);

    double v = dot(a[0], cross(a[1], a[2]));   //  volume

    for (int i = 0; i < b.size(); i++)
    {
      if (i > Nd - 1)
      {
        a[i].clear();
        b[i].clear();
      }
      else
      {
        for (int j = 0; j < b[i].size(); j++)
        {
          b[i][j] = 2. * M_PI * b[i][j] / v;
        }
      }
    }
    a.resize(Nd);
    b.resize(Nd);
  }


  void calc_element(const vector<CPPL::zcovector> &vec, const CPPL::zhematrix &mat, vector<complex<double>> &op)
  {
    op.resize(vec.size() - 1);
    {
      CPPL::zcovector tv = mat * vec[0];
      for (int j = 1; j < vec.size(); j++)
      {
        complex<double> sum(0, 0);
        for (int i = 0; i < vec[j].l; i++)
        {
          sum += conj(vec[j](i)) * tv(i);
        }
        op[j - 1] = sum;
      }
    }
  }

  // ΔS 要素の計算
  void calc_del(const vector<CPPL::zcovector> &vec, const CPPL::zhematrix &mat, vector<vector<complex<double>>> &op) 
  {
    op.resize(vec.size() - 1);
    for (int i = 0; i < op.size(); i++)
      op[i].resize(op.size());
    {
      vector<CPPL::zcovector> tv(vec.size());
      for (int i = 0; i < vec.size(); i++)
        tv[i] = mat * vec[i];

      for (int j = 1; j < vec.size(); j++)
      {
        for (int k = 1; k < vec.size(); k++)
        {
          complex<double> expec_es(0, 0);
          complex<double> expec_gs(0, 0);
          for (int i = 0; i < vec[j].l; i++)
          {
            expec_es += conj(vec[j](i)) * tv[k](i);
            expec_gs += conj(vec[0](i)) * tv[0](i);
          }

          op[j - 1][k - 1] = expec_es;
          if (j == k)
            op[j - 1][k - 1] += -expec_gs;  // δjk
        }
      }
    }
  }

  vector<vector<complex<double>>> OP;
  vector<vector<vector<complex<double>>>> delP;

  void eval_op()
  {
    OP.resize(N);
    delP.resize(N);
    for (int i = 0; i < N; ++i)
    {
      calc_element(v[i / Nc], MF[i].op(), OP[i]);
      calc_del(v[i / Nc], MF[i].op(), delP[i]);
    }
  }

public:
  mfsw(int N_, int Nc_, int Nb_, int Nitr_, int Nr_, double EPS_ = 1e-20, double MIX_ = 0.05) : N(N_), Nc(Nc_), M(N_ / Nc_), Nb(Nb_), Nr(Nr_),
                                                                                            bond(Nb_), hs(N_, 0), MF(N_), ave(N_),
                                                                                            Nitr(Nitr_), MIX(MIX_), EPS(EPS_), ene_min(1e8) {}

  tuple<int, int, exch_type> &set_J(const int n)
  {
    // exch_type>> Js is Nc x Nc matrix
    return bond[n];
  }
  void reci_v(vector<vector<double>> &a, vector<vector<double>> &b)
  { reciprocal(a, b); }
  double &set_H(const int i) { return hs[i]; }
  double &set_T() { return T; }
  CPPL::zhematrix &set_mat(const int n)
  {
    return MF[n].assign(MIX);
  }

  double mf_val(const int n) { return ave[n]; }
  double sw_ene(int i) { return ene[i]; }
  double mf_ene() { return ene_min; }
  complex<double> Vint(const int c3, const int c2, const int c1) { return Vbint[c3][c2][c1]; }
  complex<double> sw_Drate(const int i) {return Drate[i];}
  complex<double> sw_Dself(const int i) {return Dself[i];}

  string mf_out()
  {
    std::stringstream ss;
    ss << boost::format(" %12.5e   ") % mf_ene() << endl;
    for (int i = 0; i < N; i++)
    {
      //ss << ", " <<  mf_val(i);
      ss << std::fixed << std::setprecision(10) << mf_val(i) << " ";
    } 

    return ss.str();
  }

  void exec_mf(const string &fn) // 初期値をファイルから入力
  {
    vector<vector<double>> cand;
    read_file(fn, cand);

    int Ncand = cand.size();
    for (int i = 0; i < Ncand; ++i)
      exec_mf(cand[i]);
  }
  void exec_mf()  // 初期値をランダムに生成
  {
    vector<double> cand(N);

    int Ncand = cand.size();

    random_device srd;
    RND rnd(srd());
    for (int r = 0; r < Nr; r++)
    {
      cand.assign(N, 0.);
      for (int i = 0; i < Ncand; i++)
      {
          cand[i] = rnd.d(2.0)-1.0;
      }
      exec_mf(cand);
    }
  }

  void exec_mf(vector<double> &init)
  {
    for (int j = 0; j < N; j++)
    {
      MF[j].set_init(init[j]);
    }

    int Ns = MF[0].op().m;

    bool conv = 0;
    double eps;

    vector<vector<double>> wtemp(M);          
    vector<vector<CPPL::zcovector>> vtemp(M);

    for (int i = 0; i < Nitr; i++)
    {

      vector<CPPL::zhematrix> H(M);
      for (int j = 0; j < M; ++j)
      {
         wtemp[j].clear();
         vtemp[j].clear();
         H[j].resize(Ns);
         H[j].zero();
      }

      for (int j = 0; j < M; ++j)
      {
        for (int k = 0; k < Nc; ++k)
        {
          H[j] = H[j] - hs[k + Nc * j] * MF[k + Nc * j].op();
        }
      }

      for (int l = 0; l < bond.size(); ++l)
      {
        int j1 = get<0>(bond[l]);
        int j2 = get<1>(bond[l]);
        const exch_type &Js = get<2>(bond[l]);
        for (int k1 = 0; k1 < Nc; ++k1)
        {
          for (int k2 = 0; k2 < Nc; ++k2)
          {
            H[j1] = H[j1] + Js(k1, k2) * MF[k2 + Nc * j2].ave() * MF[k1 + Nc * j1].op();
            H[j2] = H[j2] + Js(k1, k2) * MF[k1 + Nc * j1].ave() * MF[k2 + Nc * j2].op();
          }
        }
      }

      for (int j = 0; j < M; ++j)
      {
        H[j].zheev(wtemp[j], vtemp[j]);
      }

      eps = 0;
      for (int j = 0; j < N; ++j)
      {
        eps += MF[j].calc(vtemp[j / Nc][0]);
      }

      if (eps < EPS)
      {
        conv = 1;
        break;
      }
    }

    double total_ene = 0;
    for (int j = 0; j < M; ++j)
      total_ene += wtemp[j][0];

    for (int l = 0; l < bond.size(); ++l)
    {
      int j1 = get<0>(bond[l]);
      int j2 = get<1>(bond[l]);
      const exch_type &Js = get<2>(bond[l]);
      for (int k1 = 0; k1 < Nc; ++k1)
      {
        for (int k2 = 0; k2 < Nc; ++k2)
        {
          total_ene += -Js(k1, k2) * MF[k2 + Nc * j2].ave() * MF[k1 + Nc * j1].ave();
        }
      }
    }
    if (conv)
    {
      double temp_ene = total_ene / double(M);
      if (temp_ene < ene_min)
      {
        ene_min = temp_ene;
        for (int i = 0; i < ave.size(); ++i)
        {
          ave[i] = MF[i].ave();
        }
        w = wtemp;
        v = vtemp;
      }
    }
  }

  /* num = 2 で k, -k のエネルギーを格納 */
  vector<CPPL::zhematrix> exec_sw(const vector<complex<double>> &gamma, const int &num = 1)
  {
    if (OP.size() == 0)
    {
      eval_op();
    }

    int Ne = OP[0].size(); // 局所状態の数-1

    int Nm = Ne * M;

    ene.assign(num * Nm, 0.);

    CPPL::zhematrix A(Nm);      A.zero();
    CPPL::zgematrix B(Nm, Nm);  B.zero();
    CPPL::zhematrix C(Nm);      C.zero();

    for (int k = 0; k < M; ++k)
    {
      for (int n = 0; n < Ne; ++n)
      {
        int j = k * Ne + n;
        A(j, j) = w[k][n + 1] - w[k][0];
        C(j, j) = w[k][n + 1] - w[k][0];
      }
    }

    for (int l = 0; l < Nb; ++l)
    {
      int j1 = get<0>(bond[l]);
      int j2 = get<1>(bond[l]);
      const exch_type &Js = get<2>(bond[l]);
      for (int n1 = 0; n1 < Ne; ++n1)
      {
        for (int n2 = 0; n2 < Ne; ++n2)
        {
          complex<double> Jtt(0), JttmT(0), Jto(0), Jot(0);
          
          for (int k1 = 0; k1 < Nc; ++k1)
          {
            for (int k2 = 0; k2 < Nc; ++k2)
            {
              Jtt += Js(k1, k2) * gamma[l] * OP[k1 + Nc * j1][n1] * OP[k2 + Nc * j2][n2];
              JttmT += Js(k1, k2) * conj(gamma[l]) * OP[k1 + Nc * j1][n2] * OP[k2 + Nc * j2][n1];
              Jto += Js(k1, k2) * gamma[l] * OP[k1 + Nc * j1][n1] * conj(OP[k2 + Nc * j2][n2]);
              Jot += Js(k1, k2) * gamma[l] * conj(OP[k1 + Nc * j1][n1]) * OP[k2 + Nc * j2][n2];
            }
          }
          A(j1 * Ne + n1, j2 * Ne + n2) += Jto;
          C(j1 * Ne + n1, j2 * Ne + n2) += Jot;
          B(j1 * Ne + n1, j2 * Ne + n2) += Jtt;
          B(j2 * Ne + n1, j1 * Ne + n2) += JttmT;
        }
      }
    }

    bogoliubov<complex<double>> SW;
    SW.calc(A, B, C);
    // CPPL::zhematrix Hsw(2 * Nm);
    // Hsw = SW.mcalc(A, B, C);

    sw_U.resize(Nm, Nm);         sw_U.zero();
    sw_V.resize(Nm, Nm);         sw_V.zero();
    sw_Jinv.resize(2*Nm, 2*Nm);  sw_Jinv.zero();
    
    for (int i = 0; i < ene.size(); i++)
      ene[i] = SW.eg(i);

    for (int i = 0; i < Nm; i++)
    {
      for (int j = 0 ;j < Nm; j++)
      {
        sw_U(i, j) = SW.Uel(i, j);
        sw_V(i, j) = SW.Vel(i, j);
      }
    }
    sw_Jinv = SW.Jinv_();

    CPPL::zgematrix spec(Nc, Nm);
    spec.zero();
    for (int i = 0; i < Nc; i++)
    {
      for (int j = 0; j < Nm; j++)
      {
        for (int n = 0; n < Ne; ++n)
        {
          for (int k = 0; k < M; ++k)
          {
            spec(i, j) += OP[i + Nc * k][n] * SW.Vel(k * Ne + n, j) + conj(OP[i + Nc * k][n]) * SW.Uel(k * Ne + n, j);
          }
        }
      }
    }

    vector<CPPL::zhematrix> W;

    W.assign(Nm, CPPL::zhematrix(Nc));
    for (int k = 0; k < Nm; k++)
    {
      for (int i = 0; i < Nc; i++)
      {
        for (int j = 0; j < i; j++)
        {
          W[k](i, j) = spec(i, k) * conj(spec(j, k)) / double(M);
        }
        W[k](i, i) = spec(i, k) * conj(spec(i, k)) / double(M);
      }
    }

    return W;
  }


  vector<vector<complex<double>>> exec_decay(const vector<double> &k, const vector<vector<double>> &xs,
                                             const vector<std::function<complex<double>(vector<double>, vector<double>)>> &g, 
                                             vector<vector<double>> &a, const double &dk, const double &eta)
  {
    if (OP.size() == 0)
    {
      eval_op();
    }

    int Ne = OP[0].size(); // 局所状態の数-1
    int Nm = Ne * M;
    int Nwv = 3;  // k, q, k-q
    int Nv = 6;  // the number of elements of Vint array 
    vector<CPPL::zgematrix> U(Nwv);  // Uk, Uq, Ukq(k-q)
    vector<CPPL::zgematrix> V(Nwv);  // Vk, Vq, Vkq(k-q)
    vector<vector<double>> Exene(Nwv);  //enek, eneq, enekq(k-q)

    for (int i = 0; i < U.size(); i++)
    {
      U[i].resize(Nm, Nm);
      U[i].zero();
      V[i].resize(Nm, Nm);
      V[i].zero();
      Exene[i].resize(Nm);
    }

    vector<complex<double>> gamma(g.size());
    for (int i = 0; i < gamma.size(); i++)
        gamma[i] = g[i](xs[i], k);

    vector<CPPL::zhematrix> W = exec_sw(gamma);  //LSWT and bogoliubov Trans.
    U[0] = sw_U;
    V[0] = sw_V;
    for (int i = 0; i < Exene[0].size(); i++)
      Exene[0][i] = ene[i];

    Drate.assign(Nm, 0.);
    Dself.assign(Nm, 0.);
    int cntx = 0;  
    int cnty = 0;
    int cntz = 0;
    vector<vector<double>> b;
    reci_v(a, b);
    vector<double> m(3);
    for (m[2] = 0.001; m[2] < 1.001; m[2] += dk)
    {
      cntz += 1;
      cnty = 0;
      cntx = 0;
      for (m[1] = 0.001; m[1] < 1.001; m[1] += dk)
      {
        cnty += 1;
        cntx = 0;
        for (m[0] = 0.001; m[0] < 1.001; m[0] += dk)
        {    
          cntx += 1;

          int Nd = b.size();

          // 波数 q のエネルギー計算
          vector<double> q(Nd);
          for (int i = 0; i < Nd; i++)
          {
            for (int j = 0; j < b.size(); j++)
            {
              q[i] += m[j] * b[j][i];
            }
          }          

          for (int i = 0; i < g.size(); i++)
            gamma[i] = g[i](xs[i], q);
          W = exec_sw(gamma);  //LSWT and bogoliubov Trans.
          U[1] = sw_U;
          V[1] = sw_V;
          for (int i = 0; i < Exene[1].size(); i++)
            Exene[1][i] = ene[i];
          

          
          // 波数 k-q のエネルギー計算
          vector<double> kq(Nd);
          for (int i = 0; i < Nd; i++)
          {
            kq[i] = k[i] - q[i];          
          }           

          for (int i = 0; i < gamma.size(); i++)
            gamma[i] = g[i](xs[i], kq);
          W = exec_sw(gamma);  //LSWT and bogoliubov Trans.
          U[2] = sw_U;
          V[2] = sw_V;
          for (int i = 0; i < Exene[2].size(); i++)
            Exene[2][i] = ene[i];

          ene.resize(Nwv * Nm);
          for (int i = 0; i < ene.size(); i += Nwv)
          {
            ene[i] = Exene[0][i / Nwv];      // k
            ene[i + 1] = Exene[1][i / Nwv];  // q
            ene[i + 2] = Exene[2][i / Nwv];  // k-q
          }
    
          // cout << "m1 = " << m1 << "   ,  " << ene << "    " << endl;

          vector<vector<vector<vector<complex<double> > > > > Vint(Nv);


          for (int i = 0; i < Vint.size(); i++)
          {
            Vint[i].resize(Nm);
            for (int c3 = 0; c3 < Vint[i].size(); c3++)
            {
              Vint[i][c3].resize(Nm);
              for (int c2 = 0; c2 < Vint[i][c3].size(); c2++)
              {
                Vint[i][c3][c2].assign(Nm, 0.);
              }
            }
          }

          for (int w = 0; w < Nwv; w++)
          {
            vector<double> K(Nd);;
            if (w == 0)
            {
              for (int i = 0; i < Nd; i++)
                K[i] = k[i];
            }
            else if(w == 1)
            {
              for (int i = 0; i < Nd; i++)
                K[i] = q[i];
            }
            else if (w == 2)
            {
              for (int i = 0; i < Nd; i++)
                K[i] = kq[i];
            }
            for (int l = 0; l < Nb; ++l)
            {
              int j1 = get<0>(bond[l]);
              int j2 = get<1>(bond[l]);
              const exch_type &Js = get<2>(bond[l]);
              for (int n1 = 0; n1 < Ne; ++n1)
              {
                for (int n2 = 0; n2 < Ne; ++n2)
                {
                  for (int n3 = 0; n3 < Ne; ++n3)
                  {
                    complex<double> Jdtt(0), Jtdt(0), Jdoo(0), Jodo(0);
                
                    for (int k1 = 0; k1 < Nc; ++k1)
                    {
                      for (int k2 = 0; k2 < Nc; ++k2)
                      {
                        Jdtt += Js(k1, k2) * g[l](xs[l], K) * delP[k1 + Nc * j1][n1][n3] * OP[k2 + Nc * j2][n2];
                        Jtdt += Js(k1, k2) * g[l](xs[l], K) * OP[k1 + Nc * j1][n1] * delP[k2 + Nc * j2][n2][n3];
                        Jdoo += Js(k1, k2) * g[l](xs[l], K) * conj(delP[k1 + Nc * j1][n1][n3]) * conj(OP[k2 + Nc * j2][n2]);
                        Jodo += Js(k1, k2) * g[l](xs[l], K) * conj(OP[k1 + Nc * j1][n1]) * conj(delP[k2 + Nc * j2][n2][n3]);
                      }
                    }
                    // cout << "w = " << w << "   : " << Jtd << "   " << Jdo << endl;

                    for (int c3 = 0; c3 < M; c3++)
                    {
                      for (int c2 = 0; c2 < M; c2++)
                      {
                        for (int c1 = 0; c1 < M; c1++)
                        {
                          if (w == 0)
                          {
                            if (c2 == c3)
                            {
                              if (c1 == j1 && c2 == j2)
                              {
                                Vint[4][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jodo);  
                                Vint[5][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jtdt;
                              }
                              else if (c1 == j2 && c2 == j1)
                              {
                                Vint[4][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += Jdtt;  
                                Vint[5][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += conj(Jdoo);
                              }
                            }
                            else if (c1 == c3)
                            {
                              if (c1 == j1 && c2 ==j2)
                              {
                                Vint[2][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jdtt;   
                                Vint[3][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jdoo);
                              }
                              else if (c1 == j2 && c2 == j1)
                              {
                                Vint[2][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += conj(Jodo);   
                                Vint[3][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += Jtdt;
                              }
                            }
                          }
                          else if (w == 1)
                          {
                            if (c2 == c3)
                            {
                              if (c1 == j1 && c2 == j2)
                              {
                                Vint[0][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jtdt;  
                                Vint[1][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jodo);
                              }
                              else if (c1 == j2 && c2 == j1)
                              {
                                Vint[0][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += conj(Jdoo);  
                                Vint[1][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += Jdtt;
                              }
                            }
                            else if (c1 == c3)
                            {
                              if (c1 == j1 && c2 ==j2)
                              {
                                Vint[4][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jdoo);   
                                Vint[5][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jdtt;
                              }
                              else if (c1 == j2 && c2 == j1)
                              {
                                Vint[4][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += Jtdt;   
                                Vint[5][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += conj(Jodo);
                              }
                            }
                          }
                          else if (w == 2)
                          {
                            if (c2 == c3)
                            {
                              if (c1 == j1 && c2 == j2)
                              {
                                Vint[2][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jtdt;  
                                Vint[3][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jodo);
                              }
                              else if (c1 == j2 && c2 == j1)
                              {
                                Vint[2][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += conj(Jdoo);  
                                Vint[3][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += Jdtt;
                              }
                            }
                            else if (c1 == c3)
                            {
                              if (c1 == j1 && c2 ==j2)
                              {
                                Vint[0][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jdoo);  
                                Vint[1][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jdtt;
                              }
                              else if (c1 == j2 && c2 == j1)
                              {
                                Vint[0][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += Jtdt;  
                                Vint[1][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += conj(Jodo);
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          // cout << endl;
          for (int w = 0; w < Vint.size(); w++)
          {
            for (int c3 = 0; c3 < Vint[w].size(); c3++)
            {
              for (int c2 = 0; c2 < Vint[w][c3].size(); c2++)
              {
                for (int c1 = 0; c1 < Vint[w][c3][c2].size(); c1++)
                {
                  if (real(Vint[w][c3][c2][c1] * conj(Vint[w][c3][c2][c1])) > 1e3)
                  {
                    cerr << "\033[31m Error \033[m : Vint is maybe injustice value." << endl;
                    cerr << Vint[w][c3][c2][c1] << endl;
                    exit(1);
                  }
                }
              }
            }
          }
          
          Vbint.resize(Nm);
          for (size_t  n3 = 0; n3 < Nm; n3++)
          {
            Vbint[n3].resize(Nm);
            for (size_t n2 = 0; n2 < Nm; n2++)
            {
              Vbint[n3][n2].assign(Nm, 0.);
            }
          }

          for (size_t n3 = 0; n3 < Nm; n3++)
          {
            for (size_t n2 = 0; n2 < Nm; n2++)
            {
              for (size_t n1 = 0; n1 < Nm; n1++)
              {
                for (size_t c3 = 0; c3 < M; c3++)
                {
                  for (size_t c2 = 0; c2 < M; c2++)
                  {
                    for (size_t c1 = 0; c1 < M; c1++)
                    {
                      for (size_t j3 = 0; j3 < Ne; j3++)
                      {
                        for (size_t j2 = 0; j2 < Ne; j2++)
                        {
                          for (size_t j1 = 0; j1 < Ne; j1++)
                          {
                            Vbint[n3][n2][n1] += Vint[0][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1] * conj(U[1](c1 * Ne + j1, n1)) * conj(U[2](c2 * Ne + j2, n2)) * U[0](c3 * Ne + j3, n3) 
                                                                              + conj(Vint[1][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1]) * conj(V[1](c1 * Ne + j1, n1)) * conj(V[2](c2 * Ne + j2, n2)) * V[0](c3 * Ne + j3, n3);
                            Vbint[n3][n2][n1] += Vint[2][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1] * conj(U[2](c1 * Ne + j1, n2)) * V[0](c2 * Ne + j2, n3) * conj(V[1](c3 * Ne + j3, n1)) 
                                                                              + conj(Vint[3][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1]) * conj(V[2](c1 * Ne + j1, n2)) * U[0](c2 * Ne + j2, n3) * conj(U[1](c3 * Ne + j3, n1));
                            Vbint[n3][n2][n1] += Vint[4][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1] * V[0](c1 * Ne + j1, n3) * conj(U[1](c2 * Ne + j2, n1)) * conj(V[2](c3 * Ne + j3, n2)) 
                                                                              + conj(Vint[5][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1]) * U[0](c1 * Ne + j1, n3) * conj(V[1](c2 * Ne + j2, n1)) * conj(U[2](c3 * Ne + j3, n2));
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }

          // cout << endl;
          double Tem = T;
          auto fB = [&Tem](double ene) {return 1. / (exp(ene / Tem) - 1.);};
          for (size_t n = 0; n < Nm; n++)
          {
            for (size_t n2 = 0; n2 < Nm; n2++)
            {
              for (size_t n1 = 0; n1 < Nm; n1++)
              {

                Drate[n] = Drate[n] + Vbint[n][n2][n1] * conj(Vbint[n][n2][n1])  * (eta / M_PI) /
                                    ( (ene[Nwv * n + 0] - ene[Nwv * n1 + 1] - ene[Nwv * n2 + 2]) * (ene[Nwv * n + 0] - ene[Nwv * n1 + 1] - ene[Nwv * n2 + 2]) + eta * eta)
                                  * ( fB(ene[Nwv * n1 + 1]) + fB(ene[Nwv * n2 + 2]) + 1. );
                if ( abs(ene[Nwv * n + 0] - ene[Nwv * n1 + 1] - ene[Nwv * n2 + 2]) > 1e-8 )
                {
                  Dself[n] = Dself[n] + Vbint[n][n2][n1] * conj(Vbint[n][n2][n1]) 
                                      / (ene[Nwv * n + 0] - ene[Nwv * n1 + 1] - ene[Nwv * n2 + 2])
                                      * ( fB(ene[Nwv * n1 + 1]) + fB(ene[Nwv * n2 + 2]) + 1. );   
                } 
                
              }
            }
          }

        }
        if (b.size() < 2)
        {
          cnty = 1;
          cntz = 1;
          break;
        }
      }
      if (b.size() < 3)
      {
        cntz = 1;
        break;
      }
    }

    int Nu = cntx * cnty * cntz;
    for (size_t n = 0; n < Nm; n++)
    {
      cout << "n = " << n << endl;
      cout << "before : ";
      cout << Drate[n] << endl;

      cout << " after : ";
      cout << (M_PI / 2.) * (1. / double(Nu)) * Drate[n] << endl;
    }
    for (int i = 0; i < Dself.size(); i++)
    {
      Dself[i] = (1. / 2.) * (1. / double(Nu)) * Dself[i];
      // cout << "n = " << i << " :   Dself = " << Dself[i] << endl;
    }

    for (int i = 0; i < Drate.size(); i++)
    {
      Drate[i] = (M_PI / 2.) * (1. / double(Nu)) * Drate[i];
      // cout << "n = " << i << " :   Drate = " << Drate[i] << endl;
    }
    vector<vector<complex<double>>> D = {Dself, Drate};
    return D; 
  }

  void exec_iDE(vector<complex<double>> &D, const vector<double> &k, const vector<vector<double>> &xs,
                const vector<std::function<complex<double>(vector<double>, vector<double>)>> &g,
                vector<vector<double>> &a, const double &dk, const int &Ni, const double &mix)  // run using random list
  {
    Drate.clear();
    Dself.clear();

    if (OP.size() == 0)
    {
      eval_op();
    }

    int Ne = OP[0].size(); // 局所状態の数-1
    int Nm = Ne * M;
    vector<double> Dtemp_imag;

    random_device srd;
    RND rnd(srd());
    for (int r = 0; r < 2; r++)
    {
      Dtemp_imag.assign(Nm, 0.);
      for (int i = 0; i < Nm; i++)
        Dtemp_imag[i] = abs(rnd.d(2.0) - 1.0);
      D = exec_iDE(k, xs, g, a, dk, Ni, mix, Dtemp_imag);
    }
  }

  vector<complex<double>> exec_iDE(const vector<double> &k, const vector<vector<double>> &xs, 
                                           const vector<std::function<complex<double>(vector<double>, vector<double>)>> &g,
                                           vector<vector<double>> &a, const double &dk, const int &Ni, const double &mix,
                                           vector<double> &Dtemp_imag)
  {
    if (OP.size() == 0)
    {
      eval_op();
    }

    double T_min = 1e-6;
    if (T < T_min || T == T_min)
      cout << "collision term : off" << endl;       
    else if (T > T_min)
      cout << "collision term : on" << endl;


    int Ne = OP[0].size(); // 局所状態の数-1
    int Nm = Ne * M;
    int Nwv = 3;                        
    int Nv = 6;                         // the number of elements of Vint array 
    vector<CPPL::zgematrix> U(Nwv);     // Uk, Uq, Ukq
    vector<CPPL::zgematrix> V(Nwv);     // Vk, Vq, Vkq
    vector<vector<double>> Exene(Nwv);  //enek, eneq, enekq

    for (int i = 0; i < U.size(); i++)
    {
      U[i].resize(Nm, Nm);
      U[i].zero();
      V[i].resize(Nm, Nm);
      V[i].zero();
      Exene[i].resize(Nm);
    }

    vector<complex<double>> gamma(g.size());
    for (int i = 0; i < gamma.size(); i++)
        gamma[i] = g[i](xs[i], k);

    vector<CPPL::zhematrix> W = exec_sw(gamma);  //LSWT and bogoliubov Trans.
    U[0] = sw_U;
    V[0] = sw_V;
    for (int i = 0; i < Exene[0].size(); i++)
      Exene[0][i] = ene[i];
    

    if (Drate.size() == 0)
      Drate.assign(Nm, 0.);

    vector<double> Drval(Drate.size());
    for (int i = 0; i < Drval.size(); i++)
    {
      if(Drval.size() != Dtemp_imag.size())
      {
        cerr << "Error : Drate have worng size( in mfsw_off_DE.hpp::exec_iDE )." << endl;
        std::exit(1);
      }
      Drval[i] = Dtemp_imag[i];
    }

    vector<vector<double>> b;     // 逆格子ベクトル
    reci_v(a, b);

    vector<vector<vector<vector<complex<double>>>>> VB;
    vector<vector<vector<vector<complex<double>>>>> WB;
    vector<vector<vector<vector<complex<double>>>>> CB;

    double Tem = T;
    auto fB = [&Tem](double ene) {return 1. / (exp(ene / Tem) - 1.);};


    vector<int> conv(Nm, 0);
    int cntx, cnty, cntz;
    cntx = 0; cnty = 0; cntz = 0;
    vector<double> m(3);
    for (m[2] = 0.001; m[2] < 1.001; m[2] += dk)
    {
      cntz += 1;
      cntx = 0;
      cnty = 0;
      for (m[1] = 0.001; m[1] < 1.001; m[1] += dk)
      { 
        cnty += 1;
        cntx = 0;
        for (m[0] = 0.001; m[0] < 1.001; m[0] += dk)
        { 
          cntx += 1;
        }
        if (b.size() < 2)
        {
          cnty = 1;
          cntz = 1;
          break;
        }
      }
      if (b.size() < 3)
      {
        cntz = 1;
        break;
      }
    }
    int Nu = cntx * cnty * cntz;
    for (int nr = 0; nr < Ni; nr++)
    {
      vector<complex<double>> Dall_imag(Nm, 0.);
      int cnt = 0;
      for (m[2] = 0.001; m[2] < 1.001; m[2] += dk)
      {
        for (m[1] = 0.001; m[1] < 1.001; m[1] += dk)
        { 
          for (m[0] = 0.001; m[0] < 1.001; m[0] += dk)
          { 
            int Nd = b.size();

            // 波数 q のエネルギー計算
            vector<double> q(Nd);
            for (int i = 0; i < Nd; i++)
            {
              for (int j = 0; j < b.size(); j++)
              {
                q[i] += m[j] * b[j][i];
              }
            }           

            for (int i = 0; i < gamma.size(); i++)
              gamma[i] = g[i](xs[i], q);
            W = exec_sw(gamma);  //LSWT and bogoliubov Trans.
            U[1] = sw_U;
            V[1] = sw_V;
            for (int i = 0; i < Exene[1].size(); i++)
              Exene[1][i] = ene[i];

            
            // 波数 kq のエネルギー計算
            vector<double> kq(Nd);
            for (int i = 0; i < Nd; i++)
            {
              kq[i] = k[i] - q[i];          
            }           

            for (int i = 0; i < gamma.size(); i++)
              gamma[i] = g[i](xs[i], kq);
            W = exec_sw(gamma);  //LSWT and bogoliubov Trans.
            U[2] = sw_U;
            V[2] = sw_V;
            for (int i = 0; i < Exene[2].size(); i++)
              Exene[2][i] = ene[i];

            ene.resize(Nwv * Nm);
            for (int i = 0; i < ene.size(); i += Nwv)
            {
              ene[i] = Exene[0][i / Nwv];      // k
              ene[i + 1] = Exene[1][i / Nwv];  // q
              ene[i + 2] = Exene[2][i / Nwv];  // kq
            }

            if (VB.size() != Nu)
            {
              vector<vector<vector<vector<complex<double> > > > > Vint(Nv);

              for (int i = 0; i < Vint.size(); i++)
              {
                Vint[i].resize(Nm);
                for (int c3 = 0; c3 < Vint[i].size(); c3++)
                {
                  Vint[i][c3].resize(Nm);
                  for (int c2 = 0; c2 < Vint[i][c3].size(); c2++)
                  {
                    Vint[i][c3][c2].assign(Nm, 0.);
                  }
                }
              }

              for (int w = 0; w < Nwv; w++)
              {
                // double Kx, Ky;
                vector<double> K(Nd);
                if (w == 0)
                {
                  for (int i = 0; i < Nd; i++)
                    K[i] = k[i];
                  // Kx = kx; Ky = ky;
                }
                else if(w == 1)
                {
                  for (int i = 0; i < Nd; i++)
                    K[i] = q[i];
                }
                else if (w == 2)
                {
                  for (int i = 0; i < Nd; i++)
                    K[i] = kq[i];
                }
                for (int l = 0; l < Nb; ++l)
                {
                  int j1 = get<0>(bond[l]);
                  int j2 = get<1>(bond[l]);
                  const exch_type &Js = get<2>(bond[l]);
                  for (int n1 = 0; n1 < Ne; ++n1)
                  {
                    for (int n2 = 0; n2 < Ne; ++n2)
                    {
                      for (int n3 = 0; n3 < Ne; ++n3)
                      {
                        complex<double> Jdtt(0), Jtdt(0), Jdoo(0), Jodo(0);
                    
                        for (int k1 = 0; k1 < Nc; ++k1)
                        {
                          for (int k2 = 0; k2 < Nc; ++k2)
                          {
                            Jdtt += Js(k1, k2) * g[l](xs[l], K) * delP[k1 + Nc * j1][n1][n3] * OP[k2 + Nc * j2][n2];
                            Jtdt += Js(k1, k2) * g[l](xs[l], K) * OP[k1 + Nc * j1][n1] * delP[k2 + Nc * j2][n2][n3];
                            Jdoo += Js(k1, k2) * g[l](xs[l], K) * conj(delP[k1 + Nc * j1][n1][n3]) * conj(OP[k2 + Nc * j2][n2]);
                            Jodo += Js(k1, k2) * g[l](xs[l], K) * conj(OP[k1 + Nc * j1][n1]) * conj(delP[k2 + Nc * j2][n2][n3]);
                          }
                        }
                        // cout << "w = " << w << "   : " << Jtd << "   " << Jdo << endl;

                        for (int c3 = 0; c3 < M; c3++)
                        {
                          for (int c2 = 0; c2 < M; c2++)
                          {
                            for (int c1 = 0; c1 < M; c1++)
                            {
                              if (w == 0)
                              {
                                if (c2 == c3)
                                {
                                  if (c1 == j1 && c2 == j2)
                                  {
                                    Vint[4][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jodo);  
                                    Vint[5][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jtdt;
                                  }
                                  else if (c1 == j2 && c2 == j1)
                                  {
                                    Vint[4][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += Jdtt;  
                                    Vint[5][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += conj(Jdoo);
                                  }
                                }
                                else if (c1 == c3)
                                {
                                  if (c1 == j1 && c2 ==j2)
                                  {
                                    Vint[2][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jdtt;   
                                    Vint[3][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jdoo);
                                  }
                                  else if (c1 == j2 && c2 == j1)
                                  {
                                    Vint[2][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += conj(Jodo);   
                                    Vint[3][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += Jtdt;
                                  }
                                }
                              }
                              else if (w == 1)
                              {
                                if (c2 == c3)
                                {
                                  if (c1 == j1 && c2 == j2)
                                  {
                                    Vint[0][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jtdt;  
                                    Vint[1][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jodo);
                                  }
                                  else if (c1 == j2 && c2 == j1)
                                  {
                                    Vint[0][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += conj(Jdoo);  
                                    Vint[1][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += Jdtt;
                                  }
                                }
                                else if (c1 == c3)
                                {
                                  if (c1 == j1 && c2 ==j2)
                                  {
                                    Vint[4][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jdoo);   
                                    Vint[5][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jdtt;
                                  }
                                  else if (c1 == j2 && c2 == j1)
                                  {
                                    Vint[4][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += Jtdt;   
                                    Vint[5][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += conj(Jodo);
                                  }
                                }
                              }
                              else if (w == 2)
                              {
                                if (c2 == c3)
                                {
                                  if (c1 == j1 && c2 == j2)
                                  {
                                    Vint[2][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jtdt;  
                                    Vint[3][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jodo);
                                  }
                                  else if (c1 == j2 && c2 == j1)
                                  {
                                    Vint[2][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += conj(Jdoo);  
                                    Vint[3][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += Jdtt;
                                  }
                                }
                                else if (c1 == c3)
                                {
                                  if (c1 == j1 && c2 ==j2)
                                  {
                                    Vint[0][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jdoo);  
                                    Vint[1][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jdtt;
                                  }
                                  else if (c1 == j2 && c2 == j1)
                                  {
                                    Vint[0][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += Jtdt;  
                                    Vint[1][c3 * Ne + n3][c2 * Ne + n1][c1 * Ne + n2] += conj(Jodo);
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }

              for (int w = 0; w < Vint.size(); w++)
              {
                for (int c3 = 0; c3 < Vint[w].size(); c3++)
                {
                  for (int c2 = 0; c2 < Vint[w][c3].size(); c2++)
                  {
                    for (int c1 = 0; c1 < Vint[w][c3][c2].size(); c1++)
                    {
                      if (real(Vint[w][c3][c2][c1] * conj(Vint[w][c3][c2][c1])) > 1e3)
                      {
                        cerr << "\033[31mError\033[m : Vint is maybe injustice value." << endl;
                        cerr << Vint[w][c3][c2][c1] << endl;
                        exit(1);
                      }
                    }
                  }
                }
              }

              VB.resize(cnt + 1);
              VB[cnt].resize(Nm);
              Vbint.resize(Nm);
              for (size_t  n3 = 0; n3 < Nm; n3++)
              {
                VB[cnt][n3].resize(Nm);
                Vbint[n3].resize(Nm);
                for (size_t n2 = 0; n2 < Nm; n2++)
                {
                  VB[cnt][n3][n2].assign(Nm, 0.);
                  Vbint[n3][n2].assign(Nm, 0.);
                }
              }

              for (size_t n3 = 0; n3 < Nm; n3++)
              {
                for (size_t n2 = 0; n2 < Nm; n2++)
                {
                  for (size_t n1 = 0; n1 < Nm; n1++)
                  {
                    for (size_t c3 = 0; c3 < M; c3++)
                    {
                      for (size_t c2 = 0; c2 < M; c2++)
                      {
                        for (size_t c1 = 0; c1 < M; c1++)
                        {
                          for (size_t j3 = 0; j3 < Ne; j3++)
                          {
                            for (size_t j2 = 0; j2 < Ne; j2++)
                            {
                              for (size_t j1 = 0; j1 < Ne; j1++)
                              {
                                Vbint[n3][n2][n1] += Vint[0][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1] * conj(U[1](c1 * Ne + j1, n1)) * conj(U[2](c2 * Ne + j2, n2)) * U[0](c3 * Ne + j3, n3) 
                                                    + conj(Vint[1][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1]) * conj(V[1](c1 * Ne + j1, n1)) * conj(V[2](c2 * Ne + j2, n2)) * V[0](c3 * Ne + j3, n3);
                                Vbint[n3][n2][n1] += Vint[2][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1] * conj(U[2](c1 * Ne + j1, n2)) * V[0](c2 * Ne + j2, n3) * conj(V[1](c3 * Ne + j3, n1)) 
                                                    + conj(Vint[3][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1]) * conj(V[2](c1 * Ne + j1, n2)) * U[0](c2 * Ne + j2, n3) * conj(U[1](c3 * Ne + j3, n1));
                                Vbint[n3][n2][n1] += Vint[4][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1] * V[0](c1 * Ne + j1, n3) * conj(U[1](c2 * Ne + j2, n1)) * conj(V[2](c3 * Ne + j3, n2)) 
                                                    + conj(Vint[5][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1]) * U[0](c1 * Ne + j1, n3) * conj(V[1](c2 * Ne + j2, n1)) * conj(U[2](c3 * Ne + j3, n2));
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }

              for (int n = 0; n < Nm; n++)
              {
                for (size_t n2 = 0; n2 < Nm; n2++)
                {
                  for (size_t n1 = 0; n1 < Nm; n1++)
                  {
                    VB[cnt][n][n2][n1] = Vbint[n][n2][n1] * conj(Vbint[n][n2][n1]);
                  }
                }
              }
              
            }
            // cout << " Nu = " << Nu << endl;
            for (int n = 0; n < Nm; n++)
            {
              if (conv[n] == 0)
              {
                for (size_t n2 = 0; n2 < Nm; n2++)
                {
                  for (size_t n1 = 0; n1 < Nm; n1++)
                  {
                    complex<double> im(0., 1.);
                    Dall_imag[n] = Dall_imag[n] - imag( VB[cnt][n][n2][n1] / 
                                    ( (ene[Nwv * n + 0] - ene[Nwv * n1 + 1] - ene[Nwv * n2 + 2]) + im * Drval[n])
                                    * ( fB(ene[Nwv * n1 + 1]) + fB(ene[Nwv * n2 + 2]) + 1. ) );
                  }
                }
              }
            }

            if (T > T_min)
            {

              for (int i = 0; i < Nd; i++)
              {
                kq[i] = k[i] + q[i];          
              }           

              for (int i = 0; i < gamma.size(); i++)
                gamma[i] = g[i](xs[i], kq);
              W = exec_sw(gamma);  //LSWT and bogoliubov Trans.
              U[2] = sw_U;
              V[2] = sw_V;
              for (int i = 0; i < Exene[2].size(); i++)
                Exene[2][i] = ene[i];

              ene.resize(Nwv * Nm);
              for (int i = 0; i < ene.size(); i += Nwv)
              {
                ene[i] = Exene[0][i / Nwv];      // k
                ene[i + 1] = Exene[1][i / Nwv];  // q
                ene[i + 2] = Exene[2][i / Nwv];  // kq
              }
        
              // cout << "m1 = " << m1 << "   ,  " << ene << "    " << endl;

              if (CB.size() != Nu)
              {
                vector<vector<vector<vector<complex<double> > > > > Vint(Nv);

                for (int i = 0; i < Vint.size(); i++)
                {
                  Vint[i].resize(Nm);
                  for (int c3 = 0; c3 < Vint[i].size(); c3++)
                  {
                    Vint[i][c3].resize(Nm);
                    for (int c2 = 0; c2 < Vint[i][c3].size(); c2++)
                    {
                      Vint[i][c3][c2].assign(Nm, 0.);
                    }
                  }
                }

                for (int w = 0; w < Nwv; w++)
                {
                  vector<double> K(Nd);
                  // double Kx, Ky;
                  if (w == 0)
                  {
                    for (int i = 0; i < Nd; i++)
                      K[i] = k[i];
                  }
                  else if(w == 1)
                  {
                    for (int i = 0; i < Nd; i++)
                      K[i] = q[i];
                  }
                  else if (w == 2)
                  {
                    for (int i = 0; i < Nd; i++)
                      K[i] = kq[i];
                  }
                  for (int l = 0; l < Nb; ++l)
                  {
                    int j1 = get<0>(bond[l]);
                    int j2 = get<1>(bond[l]);
                    const exch_type &Js = get<2>(bond[l]);
                    for (int n1 = 0; n1 < Ne; ++n1)
                    {
                      for (int n2 = 0; n2 < Ne; ++n2)
                      {
                        for (int n3 = 0; n3 < Ne; ++n3)
                        {
                          complex<double> Jdtt(0), Jtdt(0), Jdoo(0), Jodo(0);
                      
                          for (int k1 = 0; k1 < Nc; ++k1)
                          {
                            for (int k2 = 0; k2 < Nc; ++k2)
                            {
                              Jdtt += Js(k1, k2) * g[l](xs[l], K) * delP[k1 + Nc * j1][n1][n3] * OP[k2 + Nc * j2][n2];
                              Jtdt += Js(k1, k2) * g[l](xs[l], K) * OP[k1 + Nc * j1][n1] * delP[k2 + Nc * j2][n2][n3];
                              Jdoo += Js(k1, k2) * g[l](xs[l], K) * conj(delP[k1 + Nc * j1][n1][n3]) * conj(OP[k2 + Nc * j2][n2]);
                              Jodo += Js(k1, k2) * g[l](xs[l], K) * conj(OP[k1 + Nc * j1][n1]) * conj(delP[k2 + Nc * j2][n2][n3]);
                            }
                          }
                          // cout << "w = " << w << "   : " << Jtd << "   " << Jdo << endl;

                          for (int c3 = 0; c3 < M; c3++)
                          {
                            for (int c2 = 0; c2 < M; c2++)
                            {
                              for (int c1 = 0; c1 < M; c1++)
                              {
                                if (w == 0)
                                {
                                  if (c2 == c3)
                                  {
                                    if (c1 == j1 && c2 == j2)
                                    {
                                      Vint[0][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jtdt;  
                                      Vint[1][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jodo);
                                    }
                                    else if (c1 == j2 && c2 == j1)
                                    {
                                      Vint[0][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jdoo);  
                                      Vint[1][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jdtt;
                                    }
                                  }
                                  else if (c1 == c3)
                                  {
                                    if (c1 == j1 && c2 ==j2)
                                    {
                                      Vint[4][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jdoo);   
                                      Vint[5][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jdtt;
                                    }
                                    else if (c1 == j2 && c2 == j1)
                                    {
                                      Vint[4][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jtdt;   
                                      Vint[5][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jodo);
                                    }
                                  }
                                }
                                else if (w == 1)
                                {
                                  if (c2 == c3)
                                  {
                                    if (c1 == j1 && c2 == j2)
                                    {
                                      Vint[2][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jtdt;  
                                      Vint[3][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jodo);
                                    }
                                    else if (c1 == j2 && c2 == j1)
                                    {
                                      Vint[2][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jdoo);  
                                      Vint[3][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jdtt;
                                    }
                                  }
                                  else if (c1 == c3)
                                  {
                                    if (c1 == j1 && c2 ==j2)
                                    {
                                      Vint[0][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jdoo);   
                                      Vint[1][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jdtt;
                                    }
                                    else if (c1 == j2 && c2 == j1)
                                    {
                                      Vint[0][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jtdt;   
                                      Vint[1][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jodo);
                                    }
                                  }
                                }
                                else if (w == 2)
                                {
                                  if (c2 == c3)
                                  {
                                    if (c1 == j1 && c2 == j2)
                                    {
                                      Vint[4][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jodo);  
                                      Vint[5][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jtdt;
                                    }
                                    else if (c1 == j2 && c2 == j1)
                                    {
                                      Vint[4][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jdtt;  
                                      Vint[5][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jdoo);
                                    }
                                  }
                                  else if (c1 == c3)
                                  {
                                    if (c1 == j1 && c2 ==j2)
                                    {
                                      Vint[2][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jtdt;  
                                      Vint[3][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jdoo);
                                    }
                                    else if (c1 == j2 && c2 == j1)
                                    {
                                      Vint[2][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += conj(Jodo);  
                                      Vint[3][c3 * Ne + n3][c2 * Ne + n2][c1 * Ne + n1] += Jtdt;
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
                // cout << endl;

                for (int w = 0; w < Vint.size(); w++)
                {
                  for (int c3 = 0; c3 < Vint[w].size(); c3++)
                  {
                    for (int c2 = 0; c2 < Vint[w][c3].size(); c2++)
                    {
                      for (int c1 = 0; c1 < Vint[w][c3][c2].size(); c1++)
                      {
                        if (real(Vint[w][c3][c2][c1] * conj(Vint[w][c3][c2][c1])) > 1e3)
                        {
                          cerr << "\033[31mError\033[m : Vint is maybe injustice value." << endl;
                          cerr << Vint[w][c3][c2][c1] << endl;
                          exit(1);
                        }
                      }
                    }
                  }
                }
                CB.resize(cnt + 1);
                CB[cnt].resize(Nm);
                Cbint.resize(Nm);

                for (size_t n3 = 0; n3 < Nm; n3++)
                {
                  CB[cnt][n3].resize(Nm);
                  Cbint[n3].resize(Nm);
                  for (size_t n2 = 0; n2 < Nm; n2++)
                  {
                    CB[cnt][n3][n2].assign(Nm, 0.);
                    Cbint[n3][n2].assign(Nm, 0.);
                  }
                }

                for (size_t n3 = 0; n3 < Nm; n3++)
                {
                  for (size_t n2 = 0; n2 < Nm; n2++)
                  {
                    for (size_t n1 = 0; n1 < Nm; n1++)
                    {
                      for (size_t c3 = 0; c3 < M; c3++)
                      {
                        for (size_t c2 = 0; c2 < M; c2++)
                        {
                          for (size_t c1 = 0; c1 < M; c1++)
                          {
                            for (size_t j3 = 0; j3 < Ne; j3++)
                            {
                              for (size_t j2 = 0; j2 < Ne; j2++)
                              {
                                for (size_t j1 = 0; j1 < Ne; j1++)
                                {
                                  Cbint[n3][n2][n1] += Vint[0][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1] * conj(U[0](c1 * Ne + j1, n1)) * conj(U[1](c2 * Ne + j2, n2)) * U[2](c3 * Ne + j3, n3) 
                                                      + conj(Vint[1][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1]) * conj(V[0](c1 * Ne + j1, n1)) * conj(V[1](c2 * Ne + j2, n2)) * V[2](c3 * Ne + j3, n3);
                                  Cbint[n3][n2][n1] += Vint[2][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1] * conj(U[1](c1 * Ne + j1, n2)) * V[2](c2 * Ne + j2, n3) * conj(V[0](c3 * Ne + j3, n1)) 
                                                      + conj(Vint[3][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1]) * conj(V[1](c1 * Ne + j1, n2)) * U[2](c2 * Ne + j2, n3) * conj(U[0](c3 * Ne + j3, n1));
                                  Cbint[n3][n2][n1] += Vint[4][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1] * V[2](c1 * Ne + j1, n3) * conj(U[0](c2 * Ne + j2, n1)) * conj(V[1](c3 * Ne + j3, n2)) 
                                                      + conj(Vint[5][c3 * Ne + j3][c2 * Ne + j2][c1 * Ne + j1]) * U[2](c1 * Ne + j1, n3) * conj(V[0](c2 * Ne + j2, n1)) * conj(U[1](c3 * Ne + j3, n2));
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
                for (int n = 0; n < Nm; n++)
                {
                  for (size_t n2 = 0; n2 < Nm; n2++)
                  {
                    for (size_t n1 = 0; n1 < Nm; n1++)
                    {
                      CB[cnt][n][n2][n1] = Cbint[n][n2][n1] * conj(Cbint[n][n2][n1]);
                    }
                  }
                }

              }
              // cout << CB[0][0][0][0] << endl;
              // cout << endl;
              for (size_t n = 0; n < Nm; n++)
              {
                if (conv[n] == 0)
                {
                  for (size_t n2 = 0; n2 < Nm; n2++)
                  {
                    for (size_t n1 = 0; n1 < Nm; n1++)
                    {
                      complex<double> im(0., 1.);
                      Dall_imag[n] = Dall_imag[n] - imag( 2. * CB[cnt][n2][n1][n] /  
                                                ( (ene[Nwv * n + 0] + ene[Nwv * n1 + 1] - ene[Nwv * n2 + 2]) + im * Drval[n])
                                                                    * ( fB(ene[Nwv * n1 + 1]) - fB(ene[Nwv * n2 + 2]) ) );
                    }
                  }
                }
              }

            }
            // cout << endl;
            cnt += 1;
          }
          if (b.size() < 2)
            break;
        }
        if (b.size() < 3)   
          break;
      }

      for (int n = 0; n < Dall_imag.size(); n++)
      {
        Dall_imag[n] = (1. / 2.) * (1. / double(Nu)) * Dall_imag[n];
      }
      
      vector<double> eps; 
      double Eps;
      eps.assign(Nm, 0.);
      for (int n = 0; n < Dall_imag.size(); n++)
      {
        Dtemp_imag[n] = real(Dall_imag[n]);
        if (conv[n] == 0)
        {
          eps[n] += fabs(Drval[n] - Dtemp_imag[n]);
          if (eps[n] < 1e-4)
          {
            Drate[n] = Drval[n];
            conv[n] = 1;
          }
        }
      }
      // cout << Dall_imag[0] << "    " << Dtemp[0][nr] << endl;
      int tr = 1;
      for (int n = 0; n < conv.size(); n++)
      {
        tr *= conv[n];
        // if (n == 0)
        // {
        //   cout << " η = " << n << endl;
        //   cout << "imag_part :  " << Drval[n] << "    " << Dall_imag[n] << "    , conv = " << conv[n] << endl;
        //   cout << "real_part :  " << Dsval[n] << "    " << Dall_real[n] << "    , conv = " << conv[n] << endl; 
        // }   

        if (conv[n] == 0)
        {     
          Drval[n] = (1. - mix) * Drval[n] + mix * Dtemp_imag[n];
        }
      }
      if (tr == 1)
        break;
    }
    
    // vector<vector<complex<double>>> D = {Drate};
    return Drate; 
  }

  vector<double> exec_sw_mag(const vector<vector<double>> &xs, 
                     const vector<std::function<complex<double>(vector<double>, vector<double>)>> &g,
                     vector<vector<double>> a, const double &dk)
  {
    if (OP.size() == 0)
    {
      eval_op();
    }
    int Ne = OP[0].size(); // 局所状態の数-1
    int Nm = Ne * M;

    int cntx, cnty, cntz;
    cntx = 0; cnty = 0; cntz = 0;
    vector<vector<double>> b;
    reci_v(a, b);

    vector<complex<double>> sum;
    sum.assign(N, (0.,0.));
    vector<double> m(3);
    for (m[2] = 0.001; m[2] < 1.001; m[2] += dk)
    {
      cntz += 1;
      cnty = 0;
      cntx = 0;
      for (m[1] = 0.001; m[1] < 1.001; m[1] += dk)
      {
        cnty += 1;
        cntx = 0;
        for (m[0] = 0.001; m[0] < 1.001; m[0] += dk)
        {    
          cntx += 1;

          int Nd = b.size();
          // 波数 k のエネルギー計算
          vector<double> k(Nd);
          for (int i = 0; i < Nd; i++)
          {
            for (int j = 0; j < b.size(); j++)
            {
              k[i] += m[j] * b[j][i];
            }
          }          

          vector<complex<double>> gamma(g.size());
          for (int i = 0; i < g.size(); i++)
            gamma[i] = g[i](xs[i], k);
          vector<CPPL::zhematrix> W = exec_sw(gamma, 2);

          double Tem = T;
          auto fB = [&Tem](double ene) {return 1. / (exp(ene / Tem) - 1.);};


          for (int m = 0; m < M; m++)
          {
            for (int k = 0; k < Nc; k++)
            {
              for (int i = 0; i < Ne; i++)
              {
                for (int j = 0; j < Ne; j++)
                {
                  for (int n = 0; n < Nm; n++)
                  {
                    sum[k + Nc * m] += conj(sw_Jinv(i, n)) * sw_Jinv(j, n) * fB(ene[n]) * delP[k + Nc * m][i][j];
                    sum[k + Nc * m] += conj(sw_Jinv(i, n + Nm)) * sw_Jinv(j, n + Nm) * (1. + fB(ene[n + Nm])) * delP[k + Nc * m][i][j];
                    // sum[k + Nc * m] += sw_V(i, n) * conj(sw_V(j, n)) * delP[k + Nc * m][i][j];
                  }
                }
              }
            }
          }

        }

        if (b.size() < 2)
        {
          cnty = 1;
          cntz = 1;
          break;
        }
      }

      if (b.size() < 3)
      {
        cntz = 1;
        break;
      }
    }

    vector<double> moment(N);
    int Nu = cntx * cnty * cntz;
    for (size_t i = 0; i < N; i++)
    {
      moment[i] = real(sum[i]) / double(Nu) + ave[i];
    }
    
    return moment;
  }



  string exec_sw_out(const vector<complex<double>> &gamma)
  {
    vector<CPPL::zhematrix> W = exec_sw(gamma);

    std::stringstream ss;

    for (int i = 0; i < ene.size(); i++)
    {
      ss << boost::format(" %12.5e  ") % ene[i];
    }
    ss << "    ";

    for (int i = 0; i < W[0].m; ++i)
    {
      for (int k = 0; k < W.size(); k++)
      {
        ss << boost::format(" %12.5e  ") % real(W[k](i, i));
      }
      ss << "  ";
    }

    for (int i = 0; i < W[0].m; ++i)
    {
      for (int j = i + 1; j < W[0].m; ++j)
      {
        for (int k = 0; k < W.size(); k++)
        {
          ss << boost::format(" %12.5e  ") % real(W[k](i, j));
        }
        ss << "  ";
      }
    }

    for (int i = 0; i < W[0].m; ++i)
    {
      for (int j = i + 1; j < W[0].m; ++j)
      {
        for (int k = 0; k < W.size(); k++)
        {
          ss << boost::format(" %12.5e  ") % imag(W[k](i, j));
        }
        ss << "  ";
      }
    }
    return ss.str();
  }

  string exec_iDE_out(const vector<double> &k, const vector<vector<double>> xs,
                      const vector<std::function<complex<double>(vector<double>, vector<double>)>> &g,
                      vector<vector<double>> a, const double &dk, const int &Ni = 500, const double &mix = 0.5)
  {
    std::stringstream ss;
    vector<complex<double>> D;
    exec_iDE(D, k, xs, g, a, dk, Ni, mix);

    for (int i = 0; i < D.size(); i++)
    {
      ss << "   " << real(D[i]);
    }
    // cout << "   " << endl;

    return ss.str();
  }

  string exec_decay_out(const vector<double> &k, const vector<vector<double>> xs, 
                      const vector<std::function<complex<double>(vector<double>, vector<double>)>> &g,
                      vector<vector<double>> &a, const double &dk, const double &eta)
  {

    std::stringstream ss;
    vector<vector<complex<double>>> D = exec_decay(k, xs, g, a, dk, eta);

    for (int i = 0; i < D.size(); i++)
    {
      for (int j = 0; j < D[i].size(); j++)
      {
        if(abs(real(D[i][j])) > 1e3)
          ss << 0 << "   ";
        else
          ss << real(D[i][j]) << "   ";
      }
    }

      return ss.str();
  }

  string exec_sw_mag_out(const vector<vector<double>> xs, 
                      const vector<std::function<complex<double>(vector<double>, vector<double>)>> &g,
                      vector<vector<double>> &a, const double &dk)
  {
    vector<double> mag;    
    mag = exec_sw_mag(xs, g, a, dk);

    std::stringstream ss;
    for (int i = 0; i < mag.size(); i++)
    {
      ss << mag[i] << "  ";
    }
    ss << "    ";
    return ss.str();
  }



  /* ボンド自動化用 */
  void set_bond(vector<vector<double>> &xs, const int Nb, makedata<exch_type> &md)
  {
    if (M != md.make_M())
    {
      std::cerr << "\033[31mError\033[m : Can't set the bond automatically."  << endl;
      std::exit(1);
    }
    xs.resize(Nb);

    vector<tuple<vector<double>, int, int, exch_type>> Jset_sl;
    md.make_Jset_sl(Jset_sl);


    int cnt = 0;
    for (int i = 0; i < Jset_sl.size(); i++)
    {
      vector<double> delxi = std::get<0>(Jset_sl[i]);
      int mi1 = std::get<1>(Jset_sl[i]);
      int mi2 = std::get<2>(Jset_sl[i]);
      for (int j = i; j < Jset_sl.size(); j++)
      {
        vector<double> delxj = std::get<0>(Jset_sl[j]);
        int mj1 = std::get<1>(Jset_sl[j]);
        int mj2 = std::get<2>(Jset_sl[j]);
        if (j != i)
        {
          double eps = 1e-4;
          if (((mj1 == mi1) && (mj2 == mi2)) && ((mj1 != mi2) && (mj2 != mi1)))
          {
            double diff = 0.;
            for (int k = 0; k < delxi.size(); k++)
              diff = diff + abs(delxi[k] - delxj[k]);

            if(diff < eps)
            {
              double sum = 0.;
              exch_type Ji = std::get<3>(Jset_sl[i]);
              exch_type Jj = std::get<3>(Jset_sl[j]);
              for (int k1 = 0; k1 < Ji.n; k1++)
              {
                for (int k2 = 0; k2 < Ji.m; k2++)
                {
                  sum = sum + abs(Ji(k1, k2) - Jj(k1, k2));
                }
              }
              if (sum < eps)
              {
                Jset_sl.erase(Jset_sl.begin() + j);
                break;
              }
            }
          }

          if (((mj1 == mi1) && (mj2 == mi2)) && ((mj1 == mi2) && (mj2 == mi1)))  // 全てのラベルが同じ
          {
            double diff1 = 0.;
            double diff2 = 0.;
            for (int k = 0; k < delxi.size(); k++)
            {
              diff1 = diff1 + abs(delxi[k] - delxj[k]);
              diff2 = diff2 + abs(delxi[k] + delxj[k]);    // 逆位相
            }
            if(diff1 < eps)
            {
              double sum = 0.;
              exch_type Ji = std::get<3>(Jset_sl[i]);
              exch_type Jj = std::get<3>(Jset_sl[j]);
              for (int k1 = 0; k1 < Ji.n; k1++)
              {
                for (int k2 = 0; k2 < Ji.m; k2++)
                {
                  sum = sum + abs(Ji(k1, k2) - Jj(k1, k2));
                }
              }
              if (sum < eps)
              {
                Jset_sl.erase(Jset_sl.begin() + j);
                break;
              }
            }
            else if(diff2 < eps)
            {
              double sum = 0.;
              exch_type Ji = std::get<3>(Jset_sl[i]);
              exch_type Jj = std::get<3>(Jset_sl[j]);
              for (int k1 = 0; k1 < Ji.n; k1++)
              {
                for (int k2 = 0; k2 < Ji.m; k2++)
                {
                  sum = sum + abs(Ji(k1, k2) - Jj(k2, k1));     // anti-symmetric matrix まで考慮
                }
              }
              if (sum < eps)
              {
                Jset_sl.erase(Jset_sl.begin() + j);
                break;
              }
            }
          }
        }

      }
    }
    
    if (Jset_sl.size() != Nb)
    {
      std::cerr << "\033[31mError\033[m : Can't set the bond automatically."  << endl;
      std::exit(1);
    }

    // cout << "# auto bond" << endl;
    // for (int i = 0; i < Jset_sl.size(); i++)
    // {
    //   cout << std::get<0>(Jset_sl[i]) << "  ";
    //   // cout << endl;
    //   cout << std::get<1>(Jset_sl[i]) << "  ";
    //   // cout << endl;
    //   cout << std::get<2>(Jset_sl[i]) << "  ";
    //   cout << endl;
    //   cout << std::get<3>(Jset_sl[i]) << "  " << endl;
    // }

    for (int n = 0; n < Jset_sl.size(); n++)
    {
      xs[n] = std::get<0>(Jset_sl[n]);
      int i = std::get<1>(Jset_sl[n]);
      int j = std::get<2>(Jset_sl[n]);
      exch_type J = std::get<3>(Jset_sl[n]);
      set_J(n) = std::tuple<int, int, exch_type>{i, j, J};
    }
  }


};

#endif
