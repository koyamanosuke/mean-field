#ifndef MFSW_HPP
#define MFSW_HPP


#include <iostream>
#include <complex>
#include <boost/format.hpp>
#include "cpplapack/cpplapack.h"
#include <iterator>

using namespace std;

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

class mfsw
{

private:
  const int Nitr;
  const double MIX;
  const double EPS;

  int N;  // 平均場の数
  int M;  // 副格子の数
  int Nc; // 副格子あたりの平均場の数
  int Nb; //ボンドの種類

  vector<tuple<int, int, CPPL::dsymatrix>> bond; // i,jサイトを繋ぐ相互作用行列
  vector<double> hs;                             // iサイトに働く外場
  vector<mf> MF;

  vector<double> ave;                // 平均場の期待値
  vector<vector<double>> w;          // 平均場ハミルトニアンの固有値
  vector<vector<CPPL::zcovector>> v; //固有ベクトル

  CPPL::zhematrix dxA;
  CPPL::zhematrix A;

  CPPL::zgematrix Vx; //速さ行列
  CPPL::zgematrix Vy;
  CPPL::zhematrix Hsw;
  CPPL::zhematrix dMx;  //H微分
  CPPL::zhematrix dMy;  //H微分
  CPPL::zgematrix omega;

  complex<double> qh;


  vector<double> ene;
  double ene_min;

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

  vector<vector<complex<double>>> OP;

  void eval_op()
  {
    OP.resize(N);
    for (int i = 0; i < N; ++i)
    {
      calc_element(v[i / Nc], MF[i].op(), OP[i]);
    }
  }

public:
  mfsw(int N_, int Nc_, int Nb_, int Nitr_, double EPS_ = 1e-13, double MIX_ = 0.05) : N(N_), Nc(Nc_), M(N_ / Nc_), Nb(Nb_),
                                                                                            bond(Nb_), hs(N_, 0), MF(N_), ave(N_), OP(0),
                                                                                            Nitr(Nitr_), MIX(MIX_), EPS(EPS_), ene_min(1e10) {}

  tuple<int, int, CPPL::dsymatrix> &set_J(const int n)
  {
    // CPPL::dsymatrix>> Js is Nc x Nc matrix
    return bond[n];
  }
  double &set_H(const int i) { return hs[i]; }

  CPPL::zhematrix &set_mat(const int n)
  {
    return MF[n].assign(MIX);
  }

  double mf_val(const int n)
  {
    return ave[n];
  }

  double mf_ene()
  {
    return ene_min;
  }

  double mf_diff()//(1,1,1)方向の相図の計算の時のみに使う
  {
    double diff = 0;
    for(int i = 0; i<N; i++)
    {
      diff += abs(ave[i] - 0.288675);   
    }
    return diff;
  }
  double mf_ave(int i)  //aveを他のファイルで使う
  {
    return ave[i];
  }

  string mf_out()
  {

    std::stringstream ss;
    ss << boost::format(" %12.5e   ") % mf_ene() << endl;
    for (int i = 0; i < N; i++)
    {
      //ss << ", " <<  mf_val(i);
      ss << mf_val(i) << " ";
    }
    return ss.str();
  }

  string spin()
  {
      std::stringstream ss;
    for (int i = 0; i < N; i++)
    {
      //ss << ", " <<  mf_val(i);
        ss << std::fixed << setprecision(10) << mf_val(i) << " ";
    }
    return ss.str();

  }

  string diff_out()
  {
    std::stringstream difference;
    difference << mf_diff();

    return difference.str();
  }

  void exec_mf(const string &fn)
  {

    vector<vector<double>> cand;
    read_file(fn, cand);

    int Ncand = cand.size();
    for (int i = 0; i < Ncand; ++i)
      exec_mf(cand[i]);
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

    vector<vector<double>> wtemp(M);          // 平均場ハミルトニアンの固有値
    vector<vector<CPPL::zcovector>> vtemp(M); //固有ベクトル

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
        const CPPL::dsymatrix &Js = get<2>(bond[l]);
        for (int k1 = 0; k1 < Nc; ++k1)
        {
          for (int k2 = 0; k2 < Nc; ++k2)
          {
            H[j1] = H[j1] + Js(k1, k2) * MF[k1 + Nc * j2].ave() * MF[k2 + Nc * j1].op();
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
      const CPPL::dsymatrix &Js = get<2>(bond[l]);
      for (int k1 = 0; k1 < Nc; ++k1)
      {
        for (int k2 = 0; k2 < Nc; ++k2)
        {
          total_ene += -Js(k1, k2) * MF[k1 + Nc * j2].ave() * MF[k2 + Nc * j1].ave();
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
};

#endif
