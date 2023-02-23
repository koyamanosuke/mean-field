#ifndef MAKEDATA_HPP
#define MAKEDATA_HPP

//////////////////////////////

// v. 1.1
// 2023 / 2 / 16
// Change interaction matrix type to template

// v 1.0
// 2023 / 2 / 3
//////////////////////////////

#include <iostream>
#include <complex>
#include <iterator>
#include <random>

#include "cpplapack/cpplapack.h"


using namespace std;

// Batch display
template<typename T, typename _Alloc = std::allocator<T>> 
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
  for (int i = 0; i < v.size(); i++)
  os << v[i] << " ";
  return os;
}

template<typename S>
class makedata
{
private:
  std::string fn_lattice;
  std::string fn_iij;
  std::string fn_base;
  vector<vector<double>> lattice;
  vector<vector<double>> iij;
  vector<vector<double>> base;

  vector<vector<double>> a;      // 　基本並進ベクトル
  vector<vector<double>> d;      // 副格子 λ のユニットセル内での位置ベクトル

  int Ns;     // 局所状態の数  ( generator の次元)
  int Nc;     // サイトあたりの平均場の数  ( SU(Nc) )
  int M;      // データファイルの副格子の数
  int N;      // データファイルの平均場の数
  int Nd;     // 実空間の次元

  using exch_type = S;

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

public:
  makedata(std::string Fn_Lattice, std::string Fn_Iij, std::string Fn_Base) : fn_lattice(Fn_Lattice), fn_iij(Fn_Iij),
                                                                             fn_base(Fn_Base)
  {
    read_file(fn_lattice, lattice);
    read_file(fn_iij, iij);
    read_file(fn_base, base);
    Ns = base.back()[2];     
    Nc = base.back()[1];
    M = base.back()[0];
    N = M * Nc; 
    Nd =  int(sqrt(lattice[0].size())); 
  }

  int &make_Ns() {return Ns;}
  int &make_Nc() {return Nc;}
  int &make_M() {return M;}

  // void make_ad(vector<vector<double>> &as, vector<vector<double>> &ds)
  // {
  //   if (a.size() == 0 || d.size() == 0)
  //     make_unitvec();

  //   as.resize(Nd);
  //   ds.resize(M);
  //   for (int i = 0; i < Nd; i++)
  //   {
  //     as[i].resize(Nd);
  //     if (i < M)
  //       ds[i].resize(Nd);
  //     for (int j = 0; j < Nd; j++)
  //     {
  //       as[i][j] = a[i][j];
  //       if (i < M)
  //         ds[i][j] = d[i][j];
  //     }
  //   }
  // }

  void make_unitvec(vector<vector<double>> &as)
  {
    if (a.size() == 0)
      make_unitvec();

    as.resize(Nd);
    for (int i = 0; i < Nd; i++)
    {
      as[i].resize(Nd);
      for (int j = 0; j < Nd; j++)
        as[i][j] = a[i][j];
    }
  }
  
  vector<CPPL::zhematrix> make_generator()
  {
    complex<double> im(0.,1.);
    vector<CPPL::zhematrix> op(Nc);  
    for (int i = 0; i < Nc; i++)
    {
      op[i].resize(Ns);
      op[i].zero();
    }

    {
      int cnt = 0;
      for (int i = 0; i < Nc; i++)
      {
        for (int m = 0; m < Ns; m++)
        {
          for (int n = 0; n < Ns; n++)
          {
            op[i](m, n) = base[n + cnt * Ns][4] + im * base[n + cnt * Ns][5];
          }
          cnt += 1;
        }
      }
    }
    return op;
  }

  void make_unitvec()
  {
    a.resize(Nd);
    for (int i = 0; i < Nd; i++)
    {
      a[i].resize(Nd);
      for (int j = 0; j < Nd; j++)
      {
        a[i][j] = lattice[0][Nd * i + j];
      }
    }

    d.resize(M);
    for (int i = 0; i < M; i++)
    {
      d[i].resize(Nd);
      for (int j = 0; j < Nd; j++)
      {
        d[i][j] = lattice[i + 1][j];
      }
    }
  }

  vector<tuple<vector<double>, int, int, exch_type>> Jset_sl;   // 副格子まで含んだ Jset
  void make_Jset_sl(vector<tuple<vector<double>, int, int, exch_type>> &Jset_SL) 
  {
    Jset_SL.resize(Jset_sl.size());
    for (size_t i = 0; i < Jset_sl.size(); i++)
      Jset_SL[i] = Jset_sl[i];
  }

  void make_J(vector<tuple<vector<double>, exch_type>> &Jset)
  {
    if (a.size() == 0 || d.size() == 0)
    {
      make_unitvec();
    }
  
    std::vector<vector<double>> iij_tr;      // 相互作用が有限の値を持つデータのみ格納する配列
    {
      int cnt = 0;
      for (size_t i = 2; i < iij.size(); i ++)
      {
        if (int(iij[i][0]) == 0 &&
            int(iij[i][1]) == 0 &&
            int(iij[i][2]) == 0 &&
            int(iij[i][3]) == int(iij[i][4])) {}        // 磁場項の除外
        else if (abs(iij[i][7]) > 1e-4)
        {
          iij_tr.resize(cnt+1);
          iij_tr[cnt].resize(iij[i].size());
          for (int j = 0; j < iij[i].size(); j++)
            iij_tr[cnt][j] = iij[i][j];
          cnt += 1;
        }
      }
    }

    vector<exch_type> J;
    {
      int cnt = 0;
      for (int i = 0; i < iij_tr.size(); i++)
      {
        vector<double> x1(a[0].size());
        vector<double> x2(a[0].size());
        vector<double> delx(a[0].size());
        for (int j = 0; j < x1.size(); j++)
        {
          double l1 = double(iij_tr[i][3] - 1);
          double l2 = double(iij_tr[i][4] - 1);
          x1[j] = d[l1][j];
          for (int k = 0; k < a.size(); k++)
          {
            double n = double(iij_tr[i][k]);
            x2[j] += n * a[k][j];
          }
          x2[j] = x2[j] + d[l2][j];
          delx[j] = x2[j] - x1[j];
        }
        if (i == 0)
        {
          J.resize(cnt+1);    
          J[cnt].resize(Nc, Nc);   // for dgematrix 
          J[cnt].zero();
          Jset_sl.resize(cnt+1);   
          Jset.resize(cnt+1);   

          std::get<0>(Jset[cnt]) = delx;  

          if (int(iij_tr[i][3]) < int(iij_tr[i][4]))
          {
            std::get<0>(Jset_sl[cnt]) = delx;
            std::get<1>(Jset_sl[cnt]) = int(iij_tr[i][3]-1);
            std::get<2>(Jset_sl[cnt]) = int(iij_tr[i][4]-1);
          }
          else if (int(iij_tr[i][3]) > int(iij_tr[i][4]))
          {
            for (int l = 0; l < delx.size(); l++)
              delx[l] = -delx[l];
            std::get<0>(Jset_sl[cnt]) = delx;
            std::get<1>(Jset_sl[cnt]) = int(iij_tr[i][4]-1);
            std::get<2>(Jset_sl[cnt]) = int(iij_tr[i][3]-1);
          }
          else
          {
            std::get<0>(Jset_sl[cnt]) = delx;
            std::get<1>(Jset_sl[cnt]) = int(iij_tr[i][3]-1);
            std::get<2>(Jset_sl[cnt]) = int(iij_tr[i][4]-1);
          }

          cnt += 1;
        }
        else if (i > 0 &&
                (abs(iij_tr[i][0] - iij_tr[i-1][0]) != 0 ||
                abs(iij_tr[i][1] - iij_tr[i-1][1]) != 0 ||
                abs(iij_tr[i][2] - iij_tr[i-1][2]) != 0 ||
                abs(iij_tr[i][3] - iij_tr[i-1][3]) != 0 ||
                abs(iij_tr[i][4] - iij_tr[i-1][4]) != 0))
        {
          J.resize(cnt+1);    
          J[cnt].resize(Nc, Nc);  // for dgematrix  
          J[cnt].zero();
          Jset_sl.resize(cnt+1); 
          Jset.resize(cnt+1);    

          std::get<0>(Jset[cnt]) = delx;

          // if (int(iij_tr[i][3]) < int(iij_tr[i][4]))
          // {
          //   std::get<0>(Jset_sl[cnt]) = delx;
          //   std::get<1>(Jset_sl[cnt]) = int(iij_tr[i][3]-1);
          //   std::get<2>(Jset_sl[cnt]) = int(iij_tr[i][4]-1);
          // }
          // else if (int(iij_tr[i][3]) > int(iij_tr[i][4]))
          // {
          //   for (int l = 0; l < delx.size(); l++)
          //     delx[l] = -delx[l];
          //   std::get<0>(Jset_sl[cnt]) = delx;
          //   std::get<1>(Jset_sl[cnt]) = int(iij_tr[i][4]-1);
          //   std::get<2>(Jset_sl[cnt]) = int(iij_tr[i][3]-1);
          // }
          // else
          // {
            std::get<0>(Jset_sl[cnt]) = delx;
            std::get<1>(Jset_sl[cnt]) = int(iij_tr[i][3]-1);
            std::get<2>(Jset_sl[cnt]) = int(iij_tr[i][4]-1);
          // }

          cnt += 1;
        }

        int m = iij_tr[i][5] - 1;
        int n = iij_tr[i][6] - 1;
        J[cnt-1](m, n) = iij_tr[i][7];
      }
    }

    for (int i = 0; i < J.size(); i++)
    {
      std::get<3>(Jset_sl[i]) = J[i];
      std::get<1>(Jset[i]) =J[i];
    }

    for (int i = 0; i < J.size(); i++)
    {
      int j1 = std::get<1>(Jset_sl[i]);
      int j2 = std::get<2>(Jset_sl[i]);
      vector<double> delx = std::get<0>(Jset_sl[i]);
      
      if (j1 > j2)
      {
        exch_type Jnew(J[i].n, J[i].m);
        Jnew.zero();
        for (int l = 0; l < delx.size(); l++)
          delx[l] = -delx[l];
        std::get<0>(Jset_sl[i]) = delx;
        std::get<1>(Jset_sl[i]) = j2;
        std::get<2>(Jset_sl[i]) = j1;
        for (int k1 = 0; k1 < J[i].n; k1++)
        {
          for (int k2 = 0; k2 < J[i].m; k2++)
          {
            Jnew(k1, k2) = J[i](k2, k1);
          }
        }
        std::get<3>(Jset_sl[i]) = Jnew;
      }
    }

    // for (int i = 0; i < Jset_sl.size(); i++)
    // {
    //   cout << std::get<0>(Jset_sl[i]) << " ";
    //   cout << std::get<1>(Jset_sl[i]) << " ";
    //   cout << std::get<2>(Jset_sl[i]) << " ";
    //   cout << endl;
    //   cout << std::get<3>(Jset_sl[i]) << " ";
    //   cout << endl;
    // }
  }
  
  ////* 磁場のセット *////
  void make_H(vector<vector<double>> &H, const int &Mc)
  {
    H.resize(Mc);
    for (size_t i = 0; i < M; i++)
      H[i].resize(Nc);

    for (size_t k = 0; k < iij.size(); k++)
    {
      if ( (iij[k][0] == 0 && iij[k][1] == 0 && iij[k][2] == 0) && 
          (iij[k][3] == iij[k][4]) && 
          (iij[k][5] == iij[k][6]) )
      {
        int m = iij[k][3] - 1;
        int n = iij[k][5] - 1;
        H[m][n] = iij[k][7];
      }
    }
  }

  void make_H(vector<vector<double>> &H, const int &m1, const int &m2)
  {
    if (H.size() == 0)
    {
      std::cerr << "\033[31mError\033[m : '\033[1mvoid makedata::make_H(std::vector<std::vector<double> >&, const int&, const int&)\033[m'." << endl;
      std::cerr << "   Run '\033[1mvoid makedata::make_H(std::vector<std::vector<double> >&, const int&)\033[m' first."
       << endl;
      std::exit(1);
    }

    H[m2].resize(Nc);
    for (size_t i = 0; i < Nc; i++)
      H[m2][i] = H[m1][i];
  }

  // to check make_H
  void check_H(const vector<vector<double>> &H)
  {

    if (H.size() == 0)
    {
      std::cerr << "\033[31mError\033[m : Magnetic field for all sublattices is unset." << endl;  
      std::exit(1);
    }

    for (int i = 0; i < M; i++)
    {
      if (H[i].size() == 0)
      {
        std::cerr << "\033[31mError\033[m : Magnetic field for sublattice(" 
                  << i << ") is unset." << endl;
        std::exit(1);
      }
    }

  }

  // void set_bond(const int &Mc, vector<vector<double>> &xs, const int Nb, mfsw &ms)
  // {
  //   if (Mc != M)
  //   {
  //     std::cerr << "\033[31mError\033[m : Can't set the bond automatically."  << endl;
  //     std::exit(1);
  //   }
  //   xs.resize(Nb);

  //   vector<tuple<vector<double>, int, int, exch_type>> Jset_tr;   // 同一 Jset_sl を覗く
  //   int cnt = 0;
  //   for (int i = 0; i < Jset_sl.size(); i++)
  //   {
  //     vector<double> delxi = std::get<0>(Jset_sl[i]);
  //     int mi1 = std::get<1>(Jset_sl[i]);
  //     int mi2 = std::get<2>(Jset_sl[i]);
  //     for (int j = i; j < Jset_sl.size(); j++)
  //     {
  //       vector<double> delxj = std::get<0>(Jset_sl[j]);
  //       int mj1 = std::get<1>(Jset_sl[j]);
  //       int mj2 = std::get<2>(Jset_sl[j]);
  //       if (j != i)
  //       {
  //         double eps = 1e-4;
  //         if (((mj1 == mi1) && (mj2 == mi2)) && ((mj1 != mi2) && (mj2 != mi1)))
  //         {
  //           double diff = 0.;
  //           for (int k = 0; k < delxi.size(); k++)
  //             diff = diff + abs(delxi[k] - delxj[k]);

  //           if(diff < eps)
  //           {
  //             Jset_sl.erase(Jset_sl.begin() + j);
  //             break;
  //           }
  //         }

  //         if (((mj1 == mi1) && (mj2 == mi2)) && ((mj1 == mi2) && (mj2 == mi1)))
  //         {
  //           double diff1 = 0.;
  //           double diff2 = 0.;
  //           for (int k = 0; k < delxi.size(); k++)
  //           {
  //             diff1 = diff1 + abs(delxi[k] - delxj[k]);
  //             diff2 = diff1 + abs(delxi[k] + delxj[k]);
  //           }
  //           if((diff1 < eps) || (diff2 < eps))
  //           {
  //             Jset_sl.erase(Jset_sl.begin() + j);
  //             break;
  //           }
  //         }
  //       }

  //     }
  //   }
    
  //   if (Jset_sl.size() != Nb)
  //   {
  //     std::cerr << "\033[31mError\033[m : Can't set the bond automatically."  << endl;
  //     std::exit(1);
  //   }

  //   for (int i = 0; i < Jset_sl.size(); i++)
  //   {
  //     cout << std::get<0>(Jset_sl[i]) << "  ";
  //     // cout << endl;
  //     cout << std::get<1>(Jset_sl[i]) << "  ";
  //     // cout << endl;
  //     cout << std::get<2>(Jset_sl[i]) << "  ";
  //     cout << endl;
  //     cout << std::get<3>(Jset_sl[i]) << "  " << endl;
  //   }

  //   for (int n = 0; n < Jset_sl.size(); n++)
  //   {
  //     xs[n] = std::get<0>(Jset_sl[n]);
  //     int i = std::get<1>(Jset_sl[n]);
  //     int j = std::get<2>(Jset_sl[n]);
  //     exch_type J = std::get<3>(Jset_sl[n]);
  //     ms.set_J(n) = std::tuple<int, int, exch_type>{i, j, J};
  //   }
  // }
};

#endif
