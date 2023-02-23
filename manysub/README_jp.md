# mfsw_test

`mfsw_test` は, 一般的な相互作用模型に平均場理論およびスピン波理論を適用して作成されたプログラムである `test.cpp`
を実行するテストディレクトリです。

$$\mathcal{H} \ = \ \sum_{\langle ij \rangle} \sum_{\xi\xi^\prime} I_{ij}^{\xi\xi^\prime}\mathcal{O}_i^{\xi}\mathcal{O}_j^{\xi^\prime}-\sum_{i}\sum_{\xi} H_i^\xi \mathcal{O}_{i}^\xi$$

以下を計算することができます。
- 局所状態
- スピン波分散
- 動的構造因子 $S_{\xi\xi^\prime}(\mathbf{q},\omega)$
- モーメントの縮み (or SW モーメント)
- 準粒子減衰率

動的構造因子 $S_{\xi\xi^\prime}(\mathbf{q},\omega)$ は、線形スピン波理論の範囲で次のように表されます。
$$
S_{\xi\xi^\prime} (\mathbf{q},\omega) \ = \ \sum_{\eta}^{N} \tilde{W}_{\eta,\mathbf{q}}^{\xi} \tilde{W}_{\eta,\mathbf{q}}^{\xi^\prime*}\delta (\omega - \varepsilon_{\eta,\mathbf{q}}),
$$
ここで、 $N$ はスピン波バンドの数であり 
$\varepsilon_{\eta,\mathbf{q}}$ は $\eta$ 番目のバンドの励起エネルギーです。

プログラム上ではすべての $\tilde{W}_{\eta,\mathbf{q}}^{\xi}$ が計算され、
`std::ofstream` によってオープンされたファイルに書き込まれます。

# Overview
以下は、データ `"cubic2subbcc"` から
平均場解、スピン波分散、SWモーメント、動的構造因子を計算するコードです。

```cpp:test.cpp

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

  ms.exec_mf();                    // Nr 回のランダムな初期値で平均場解を求める
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
      sw << x << " " << ms.exec_sw_out(gamma);   // スピン波エネルギーおよび動的構造因子を計算して書き込む
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
```

以下のことに注意して下さい：
### 1. ``cout << ms.mf_out() << endl;`` の出力は、次の順序で行われます：
$$
\begin{align}
&E_{\text{gs}}, \nonumber\\
&\langle\mathcal{O}^{\xi=0}\rangle_{A}, \ \langle\mathcal{O}^{\xi=1}\rangle_{A}, \cdots, \ \langle\mathcal{O}^{\xi=\xi_{\text{max}}}\rangle_{A},\ 
\langle\mathcal{O}^{\xi=0}\rangle_{B}, \ \cdots, \ \langle\mathcal{O}^{\xi=\xi_{\text{max}}}\rangle_{B}, \
\langle\mathcal{O}^{\xi=0}\rangle_{C}, \ \cdots \nonumber
\end{align}
$$
ここで、大文字のアルファベットはサブラティスを表しています。
例として、`"cubic2subbcc"`のデータに対して、以下のような出力が得られる。

```bash:output
#  -1.00000e+00   
0.7071067973 0.5845336644 -0.0056875930 -0.3978543406 0.7071066436 -0.5845336470 0.0056875929 0.3978543288 
```

### 2. ``sw << x << " " << ms.exec_sw_out(gamma);`` は以下の順序で `fnsw` に書き出します：
$$
\begin{align}
&x,\
\varepsilon_{0,\mathbf{q}},\ \varepsilon_{1,\mathbf{q}}, \cdots, \varepsilon_{N,\mathbf{q}},\ 
\Big|\tilde{W}_{0,\mathbf{q}}^{\xi=0}\Big|^2, \Big|\tilde{W}_{1,\mathbf{q}}^{\xi=0}\Big|^2 ,\cdots ,\Big|\tilde{W}_{N,\mathbf{q}}^{\xi=0}\Big|^2,\ 
\Big|\tilde{W}_{0,\mathbf{q}}^{\xi=1}\Big|^2, \Big|\tilde{W}_{1,\mathbf{q}}^{\xi=1}\Big|^2 ,\cdots ,\Big|\tilde{W}_{N,\mathbf{q}}^{\xi=1}\Big|^2,\ \cdots\cdots,
\Big|\tilde{W}_{0,\mathbf{q}}^{\xi=\xi_{\text{max}}}\Big|^2, \Big|\tilde{W}_{1,\mathbf{q}}^{\xi=\xi_{\text{max}}}\Big|^2 ,\cdots ,\Big|\tilde{W}_{N,\mathbf{q}}^{\xi=\xi_{\text{max}}}\Big|^2,\ 
\text{Re} \big[ (\text{others})\big],\ \text{Im} \big[ (\text{others})\big]
\nonumber
\end{align}
$$

### 3. ``mg << " " << ms.exec_sw_mag_out(xs, g, a, dk) << endl;`` は以下の順序で  `fnmg` に書き出します:
$$
\newcommand{\bbraket}[1]{\langle\langle#1\rangle\rangle}
\begin{align}
&\bbraket{\mathcal{O}^{\xi=0}}_{A}, \ \bbraket{\mathcal{O}^{\xi=1}}_{A}, \cdots, \ \bbraket{\mathcal{O}^{\xi=\xi_{\text{max}}}}_{A},\ 
\bbraket{\mathcal{O}^{\xi=0}}_{B}, \ \cdots, \ \bbraket{\mathcal{O}^{\xi=\xi_{\text{max}}}}_{B}, \
\bbraket{\mathcal{O}^{\xi=0}}_{C}, \ \cdots \nonumber
\end{align}
$$
ここで、 $\langle\langle ~\cdot~ \rangle\rangle$ は Bogoliubov 真空における期待値です。

# Interfaces
`makedata.hpp` および `mfsw.hpp` はそれぞれ１つのインターフェースを提供します。

### ・makedata
`makedata` クラスの主な役割は、元データから $\mathcal{O}_i^{\gamma}$, $I_{ij}$, $H_{i}$, および基本並進ベクトルを生成することです。

引数には元データファイルです。
 `"***_lattice.dat"`, `"***_iij.dat"`, `"***_base.dat"` の順に入力して下さい。
### ・mfsw
`mfsw' クラスは、計算を実行するクラスです.
平均場解では、初期値の振り方をファイル読み込みに変更することができます。
例として、以下の初期値入力ファイルを用意します。
```txt:data/input/init.txt
0.7071067812 0 0 0.7071067812  0.7071067812 0 0 -0.7071067812
```
そして `ms.exec_mf()` を次のように変更します。
```cpp:test.cpp
ms.exec_mf("data/input/init.txt");
```
プログラムを実行すると以下の出力を得ます。
``` bash:output
#  -1.00000e+00   
0.7071067812 0.0000000000 0.0000000000 0.7071067812 0.7071067812 0.0000000000 0.0000000000 -0.7071067812 
```
なお、いくつかの初期値を入れておくと、初期値の数だけ自己無撞着な方程式を解き、
その初期値の範囲内で最安定解を見つけてくれます。
# Graphics
スピン波分散と動的構造因子は `sw.py` と `spec.py` によってプロットされます。
これらの python コードはコマンドライン引数を使用しています。
```bash
$ python sw.py M Ne
$ python spec.py M Ne
```
ここで、Mは副格子の総数、Neは局所励起状態の総数です。

例として `"cubic2subbcc"` では副格子数は２、局所励起状態の数は１なので、次のように入力します：
```bash
$ python sw.py 2 1
$ python spec.py 2 1
```

このテストプログラムでは、`spec.py` は
$$
S(\mathbf{q},\omega) \ = \ \sum_{\xi} S_{\xi\xi} (\mathbf{q},\omega)
$$
のcolormapを生成します。
# Requirement
- C++11 compatible environment

- cpplapack-2015.05.11 (header-only library and present in the directory)

- boost (header-only library)