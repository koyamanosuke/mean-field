#ifndef BOGOLIUBOV_HPP
#define BOGOLIUBOV_HPP

//////////////////////////////
// v 1.1
// 2019/4/22
//
// v 1.0
// 2019/4/21
//////////////////////////////


#include "cpplapack/cpplapack.h"

class bogoliubovd{
private:
  std::size_t ndim;
  int dpotrf(CPPL::dsymatrix a, CPPL::dgematrix &b);

  std::vector<double> ev;
  CPPL::dgematrix V;
  CPPL::dgematrix U;
  CPPL::dsymatrix Am,Bm;


public:
  bogoliubovd(const std::size_t ndim_): ndim(ndim_), Am(ndim_), Bm(ndim_){}
  bogoliubovd(){}
  double eg(std::size_t i){return ev[i];}
  double to(const std::size_t i,const std::size_t j,const std::size_t k)const { return U(i,k)*U(j,k);}
  double ot(const std::size_t i,const std::size_t j, const std::size_t k)const { return V(i,k)*V(j,k);}
  double tt(const std::size_t i,std::size_t j, std::size_t k)const { return U(i,k)*V(j,k);}
  double oo(const std::size_t i,std::size_t j, std::size_t k)const { return V(i,k)*U(j,k);}
  double to(std::size_t i, std::size_t j);
  double ot(std::size_t i, std::size_t j);
  double tt(std::size_t i, std::size_t j);
  double oo(std::size_t i, std::size_t j);
  double to(std::size_t k);
  double ot(std::size_t k);
  double tt(std::size_t k);
  double oo(std::size_t k);

  std::size_t size(){return ndim;};

  CPPL::dsymatrix & getA(){return Am;};
  CPPL::dsymatrix & getB(){return Bm;};

  double & getA(const std::size_t i, const std::size_t j){return Am(i,j);};
  double & getB(const std::size_t i, const std::size_t j){return Bm(i,j);};

  int calc(const CPPL::dsymatrix A, const CPPL::dsymatrix B);
  int calc(){return calc(Am, Bm);}   //  0: success,   1: failure

};


extern "C" {
  void dpotrf_(const char &uplo, const int &N, double *A, const int &lda, int &info);
}

// the Cholesky factorization of a real symmetric positive definite matrix A
// A = L*trans(L) or trans(U)*U
// INFO    (output) INTEGER
//                = 0:  successful exit
//                < 0:  if INFO = -i, the i-th argument had an illegal value
//                > 0:  if INFO = i, the leading minor of order i is not positive
//                definite, and the factorization could not be completed.
int bogoliubovd::dpotrf(CPPL::dsymatrix a, CPPL::dgematrix &b){
  const char uplo = 'U';
  const int N = a.m;
  const int lda = a.m;
  int info;

  b.resize(N,N);
  b.zero();
  b += a;
  dpotrf_(uplo, N, b.array, lda, info);
  for(int i = 0;i<N;i++){
    for(int j = 0;j<i;j++){
      b(i,j) = 0;
    }
  }
  return info;
}

int bogoliubovd::calc(CPPL::dsymatrix A, CPPL::dsymatrix B){

  ndim = A.m;

  CPPL::dgematrix K;
  int info = dpotrf(A-B,K); //Cholesky deconposition

  CPPL::dgematrix tK = CPPL::t(K);

  CPPL::dsymatrix ApB = A+B;
  CPPL::dsymatrix C(ndim);
  {
    CPPL::dgematrix temp = K*ApB*tK;
    for(int i=0;i<ndim;i++){
      for(int j=0;j<=i;j++){
        C(i,j) = temp(i,j);
      }
    }
  }

  std::vector<double> w;
  std::vector<CPPL::dcovector> v;

  C.dsyev(w, v);

  std::vector<CPPL::dcovector> P,M;

  ev.clear();

  V.resize(ndim,ndim);
  U.resize(ndim,ndim);

  for(int i=0;i<ndim;i++){

    ev.push_back(sqrt(w[i]));  //error check

    v[i] = v[i]/sqrt(ev[i]);

    P.push_back( tK * v[i] );
    M.push_back( ApB * P[i] / ev[i] );

    for(int j=0;j<ndim;j++){
      V(j,i) = 0.5*(P[i](j) + M[i](j));
      U(j,i) = 0.5*(P[i](j) - M[i](j));
    }
  }

  return info;
}



inline double bogoliubovd::to(std::size_t i,std::size_t j){
  double sum = 0;
  for(std::size_t k=0;k<ndim;k++) sum += to(i,j,k);
  return sum;
}
inline double bogoliubovd::ot(std::size_t i,std::size_t j){
  double sum = 0;
  for(std::size_t k=0;k<ndim;k++) sum += ot(i,j,k);
  return sum;
}
inline double bogoliubovd::tt(std::size_t i,std::size_t j){
  double sum = 0;
  for(std::size_t k=0;k<ndim;k++) sum += tt(i,j,k);
  return sum;
}
inline double bogoliubovd::oo(std::size_t i,std::size_t j){
  double sum = 0;
  for(std::size_t k=0;k<ndim;k++) sum += oo(i,j,k);
  return sum;
}

inline double bogoliubovd::to(std::size_t k){
  double sum = 0;
  for(std::size_t i=0;i<ndim;i++){
    for(std::size_t j=0;j<ndim;j++){
      sum += to(i,j,k);
    }
  }
  return sum;
}
inline double bogoliubovd::ot(std::size_t k){
  double sum = 0;
  for(std::size_t i=0;i<ndim;i++){
    for(std::size_t j=0;j<ndim;j++){
      sum += ot(i,j,k);
    }
  }
  return sum;
}
inline double bogoliubovd::tt(std::size_t k){
  double sum = 0;
  for(std::size_t i=0;i<ndim;i++){
    for(std::size_t j=0;j<ndim;j++){
      sum += tt(i,j,k);
    }
  }
  return sum;
}
inline double bogoliubovd::oo(std::size_t k){
  double sum = 0;
  for(std::size_t i=0;i<ndim;i++){
    for(std::size_t j=0;j<ndim;j++){
      sum += oo(i,j,k);
    }
  }
  return sum;
}


class bogoliubovc{
private:
  std::size_t ndim;
  int zpotrf(const CPPL::zhematrix &a, CPPL::zgematrix &b);
  std::vector<double> ev;
  CPPL::zgematrix U;
  CPPL::zhematrix Am,Cm;
  CPPL::zgematrix Bm;

public:
  bogoliubovc(const std::size_t ndim_): ndim(ndim_), Am(ndim_), Bm(ndim_,ndim_), Cm(ndim_){}
  bogoliubovc(){}
  int calc(const CPPL::zhematrix& A);
  int calc(const CPPL::zhematrix& A, const CPPL::zgematrix& B, const CPPL::zhematrix& C);
  CPPL::zhematrix mcalc(const CPPL::zhematrix& A, const CPPL::zgematrix& B, const CPPL::zhematrix& C);
  double eg(std::size_t i){ return ev[i];}
  std::size_t size(){return ndim;}
  std::size_t sizef(){return ndim*2;}

  std::complex<double> Uel(std::size_t i, std::size_t j){return U(i,j);}
  std::complex<double> Vel(std::size_t i, std::size_t j){return U(i+ndim,j);}
  std::complex<double> Wel(std::size_t i, std::size_t j){return U(i,j+ndim);}
  std::complex<double> Xel(std::size_t i, std::size_t j){return U(i+ndim,j+ndim);}
  std::complex<double> Jinv(std::size_t i, std::size_t j){return U(i,j);}
  CPPL::zgematrix Jinv_(){return U;}


  std::complex<double> to(std::size_t i,std::size_t j, std::size_t k){
    return Vel(i,k)*conj(Vel(j,k));
  }
  std::complex<double> ot(std::size_t i,std::size_t j, std::size_t k){
    return Uel(i,k)*conj(Uel(j,k));
  }
  std::complex<double> tt(std::size_t i,std::size_t j, std::size_t k){
    return Vel(i,k)*conj(Uel(j,k));
  }
  std::complex<double> oo(std::size_t i,std::size_t j, std::size_t k){
    return Uel(i,k)*conj(Vel(j,k));
  }
  std::complex<double> to(std::size_t i, std::size_t j);
  std::complex<double> ot(std::size_t i, std::size_t j);
  std::complex<double> tt(std::size_t i, std::size_t j);
  std::complex<double> oo(std::size_t i, std::size_t j);
  std::complex<double> to(std::size_t k);
  std::complex<double> ot(std::size_t k);
  std::complex<double> tt(std::size_t k);
  std::complex<double> oo(std::size_t k);

  // std::complex<double> & getA(const std::size_t i, const std::size_t j){return Am(i,j);};
  std::complex<double> & getB(const std::size_t i, const std::size_t j){return Bm(i,j);};
  // std::complex<double> & getC(const std::size_t i, const std::size_t j){return Cm(i,j);};

  CPPL::zhecomplex getA(const std::size_t i, const std::size_t j){return Am(i,j);};
  CPPL::zhecomplex getC(const std::size_t i, const std::size_t j){return Cm(i,j);};

  int calc(){return calc(Am, Bm, Cm);}   //  0: success,   1: failure

};


extern "C" {
  void zpotrf_(const char &uplo, const int &N, std::complex<double> *A, const int &lda, int &info);
}
// the Cholesky factorization of a Hermite positive definite matrix A
// A = L*conjt(L) or conjt(U)*U
// INFO    (output) INTEGER
//                = 0:  successful exit
//                < 0:  if INFO = -i, the i-th argument had an illegal value
//                > 0:  if INFO = i, the leading minor of order i is not positive
//                definite, and the factorization could not be completed.


int bogoliubovc::zpotrf(const CPPL::zhematrix &a, CPPL::zgematrix &b){
  const char uplo = 'U';
  const int N = a.m;
  const int lda = a.m;
  int info;
  b.resize(N,N);
  b.zero();
  b += a;

  zpotrf_(uplo, N, b.array, lda, info);
  for(int i = 0;i<N;i++){
    for(int j = 0;j<i;j++){
      b(i,j) = 0;
    }
  }
  return info;
}


int bogoliubovc::calc(const CPPL::zhematrix& A){  
  ndim = A.m/2;
  std::size_t N = ndim;

  CPPL::zgematrix K;
  std::size_t info = zpotrf(A, K);

  CPPL::zgematrix Ip(ndim*2,ndim*2);
  Ip.zero();
  for(std::size_t i=0;i<N;i++) Ip(i,i)=1;
  for(std::size_t i=N;i<ndim*2;i++) Ip(i,i)=-1;

  CPPL::zhematrix Ch(ndim*2);
  {
    CPPL::zgematrix C = K*Ip*conjt(K);
    for(std::size_t i=0;i<ndim*2;i++){
      for(std::size_t j=i;j<ndim*2;j++){
        Ch(i,j)=C(i,j);
      }
    }
  }

  std::vector<double> wt;
  std::vector<CPPL::zcovector> vt;
  Ch.zheev(wt,vt);

  //sort
  ev.resize(ndim*2);
  std::vector<CPPL::zcovector> v = vt;
  CPPL::zgematrix Esq(ndim*2,ndim*2);
  Esq.zero();
  for(std::size_t i=0;i<N;i++){
    ev[i] = wt[i+N];
    v[i] = vt[i+N];
    ev[i+N] = -wt[N-1-i];
    v[i+N] = vt[N-1-i];
    Esq(i,i) = sqrt(ev[i]);
    Esq(i+N,i+N) = sqrt(ev[i+N]);
  }
  CPPL::zgematrix V(ndim*2,ndim*2);
  for(std::size_t i=0;i<ndim*2;i++){
    for(std::size_t j=0;j<ndim*2;j++){
      V(i,j) = v[j](i);
    }
  }
  U = i(K)*V*Esq;
  return info;
}

int bogoliubovc::calc(const CPPL::zhematrix& A, const CPPL::zgematrix& B, const CPPL::zhematrix& C){
  std::size_t N = A.m;
  CPPL::zhematrix H(2*N);
  for(std::size_t i=0;i<N;i++){
    for(std::size_t j=i;j<N;j++){
      H(i,j) = A(i,j);
      H(i+N,j+N) = C(i,j);
    }
    for(std::size_t j=0;j<N;j++){
      H(i,j+N) = B(i,j);
    }
  }
  return calc(H);
}

CPPL::zhematrix bogoliubovc::mcalc(const CPPL::zhematrix& A, const CPPL::zgematrix& B, const CPPL::zhematrix& C){
  std::size_t N = A.m;
  CPPL::zhematrix H(2*N);
  for(std::size_t i=0;i<N;i++){
    for(std::size_t j=i;j<N;j++){
      H(i,j) = A(i,j);
      H(i+N,j+N) = C(i,j);
    }
    for(std::size_t j=0;j<N;j++){
      H(i,j+N) = B(i,j);
    }
  }
  return H;
}


inline std::complex<double> bogoliubovc::to(std::size_t i,std::size_t j){
  std::complex<double> sum = 0;
  for(std::size_t k=0;k<ndim;k++) sum += to(i,j,k);
  return sum;
}
inline std::complex<double> bogoliubovc::ot(std::size_t i,std::size_t j){
  std::complex<double> sum = 0;
  for(std::size_t k=0;k<ndim;k++) sum += ot(i,j,k);
  return sum;
}
inline std::complex<double> bogoliubovc::tt(std::size_t i,std::size_t j){
  std::complex<double> sum = 0;
  for(std::size_t k=0;k<ndim;k++) sum += tt(i,j,k);
  return sum;
}
inline std::complex<double> bogoliubovc::oo(std::size_t i,std::size_t j){
  std::complex<double> sum = 0;
  for(std::size_t k=0;k<ndim;k++) sum += oo(i,j,k);
  return sum;
}

inline std::complex<double> bogoliubovc::to(std::size_t k){
  std::complex<double> sum = 0;
  for(std::size_t i=0;i<ndim;i++){
    for(std::size_t j=0;j<ndim;j++){
      sum += to(i,j,k);
    }
  }
  return sum;
}
inline std::complex<double> bogoliubovc::ot(std::size_t k){
  std::complex<double> sum = 0;
  for(std::size_t i=0;i<ndim;i++){
    for(std::size_t j=0;j<ndim;j++){
      sum += ot(i,j,k);
    }
  }
  return sum;
}
inline std::complex<double> bogoliubovc::tt(std::size_t k){
  std::complex<double> sum = 0;
  for(std::size_t i=0;i<ndim;i++){
    for(std::size_t j=0;j<ndim;j++){
      sum += tt(i,j,k);
    }
  }
  return sum;
}
inline std::complex<double> bogoliubovc::oo(std::size_t k){
  std::complex<double> sum = 0;
  for(std::size_t i=0;i<ndim;i++){
    for(std::size_t j=0;j<ndim;j++){
      sum += oo(i,j,k);
    }
  }
  return sum;
}

///////////////////////////////////////////////////////////////////////////////////
//for alias

template <typename>
struct bogoliubov_temp {
  using type = bogoliubovc;
};
template <>
struct bogoliubov_temp<double> {
  using type = bogoliubovd;
};
template <typename T>
using bogoliubov = typename bogoliubov_temp<T>::type;

///////////////////////////////////////////////////////////////////////////////////

#endif
