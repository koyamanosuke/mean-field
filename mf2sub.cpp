#include <iostream>
#include <random>
#include <complex>
#include <vector>

#include <boost/timer.hpp>
#include <boost/format.hpp>

#include "cpplapack/cpplapack.h"
#define Nitr 30000
#define mix 0.05


  using namespace std;


class mf: CPPL::zhematrix {
private:
  double val;
  CPPL::zhematrix mat;
  double init;
public:
  mf(CPPL::zhematrix mat_): mat(mat_), init(1000){val=init;}
  mf(CPPL::zhematrix mat_, double init_): mat(mat_), init(init_){val=init;}
  CPPL::zhematrix op(){
    return mat;
  }
  double calc(const CPPL::zcovector& v){
    CPPL::zcovector vt = mat*v;
    complex<double> sum(0,0);
    for(int i=0;i<v.l;i++)
      sum += conj(v(i))*vt(i);

    double val_old = val;
    val = std::real(sum)*mix + val_old*(1-mix);
    return fabs(val-val_old);
  }
  double ave(){return val;};
};

struct param{
  double K;
  double J;
  double hx;
  double hy;
  double hz;
  double Tmin;
  double Tmax;
};
void meanfield(double& ene_min, const param& pr,
               vector<double>& init, vector<double>& ave,
               const CPPL::zhematrix& Sxop, const CPPL::zhematrix& Syop, const CPPL::zhematrix& Szop,
               vector<double>& wA, vector<double>& wB,
               vector<CPPL::zcovector>& vA, vector<CPPL::zcovector>& vB)
{

  mf SxA(Sxop,init[0]);
  mf SxB(Sxop,init[1]);
  mf SyA(Syop,init[2]);
  mf SyB(Syop,init[3]);
  mf SzA(Szop,init[4]);
  mf SzB(Szop,init[5]);


  bool conv=0;
  double eps;

  double eneA, eneB;
  for(int i=0;i<Nitr;i++)
    {
      // vector<double> wA;
      // vector<double> wB;
      // vector<CPPL::zcovector> vA;
      // vector<CPPL::zcovector> vB;
      wA.clear();
      wB.clear();
      vA.clear();
      vB.clear();

    // /* Neel */
    CPPL::zhematrix HA = (2*pr.K+3*pr.J)*(SxB.ave()*SxA.op()+SyB.ave()*SyA.op()+SzB.ave()*SzA.op())
                          -pr.hx*SxA.op() - pr.hy*SyA.op() - pr.hz*SzA.op();

    CPPL::zhematrix HB = (2*pr.K+3*pr.J)*(SxA.ave()*SxB.op()+SyA.ave()*SyB.op()+SzA.ave()*SzB.op())              //Neel
                          -pr.hx*SxB.op() - pr.hy*SyB.op() - pr.hz*SzB.op();



      HA.zheev(wA,vA);
      HB.zheev(wB,vB);
      eneA = wA[0];
      eneB = wB[0];

      eps = 0;
      eps += SxA.calc(vA[0]);
      eps += SxB.calc(vB[0]);
      eps += SyA.calc(vA[0]);
      eps += SyB.calc(vB[0]);
      eps += SzA.calc(vA[0]);
      eps += SzB.calc(vB[0]);


      if(eps < 1e-8) {
          conv = 1;
          break;}
    }


  double total_ene = ( eneA + eneB - (2*pr.K+3*pr.J) * (SxA.ave() * SxB.ave()+ SyA.ave() * SyB.ave()+ SzA.ave() * SzB.ave()) 
                      - (2*pr.K+3*pr.J) * (SxB.ave() * SxA.ave()+ SyB.ave() * SyA.ave()+ SzB.ave() * SzA.ave()))/2.;   //Neel



  if(conv){
    if(total_ene < ene_min){
      ene_min = total_ene;
      ave[0] = SxA.ave();
      ave[1] = SxB.ave();
      ave[2] = SyA.ave();
      ave[3] = SyB.ave();
      ave[4] = SzA.ave();
      ave[5] = SzB.ave();

    }
  }

  // cerr << "conv = " << conv << endl;
  // cerr << "eps  = " << eps << endl;
  // cerr << "total energy  = " << total_ene << endl << endl;


}

void calc_element(const vector<CPPL::zcovector>& vec, const CPPL::zhematrix& mat, vector<complex<double> >& op)
{
  op.resize(vec.size()-1);
  {
    CPPL::zcovector tv = mat*vec[0];
    for(int j=1;j<vec.size();j++){
      complex<double> sum(0,0);
      for(int i=0;i<vec[j].l;i++){
        sum+= conj(vec[j](i))*tv(i);
      }
      op[j-1]=sum;
    }
  }
}


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




int main(int argc, char *argv[])
{

  ofstream ofs("spec2sub.txt");
  ofstream ofs1("ch2sub.txt");
  ofstream ofs2("qh2sub.txt");




  complex<double> im(0,1);

  CPPL::zhematrix Sxop(2); Sxop.zero();
  CPPL::zhematrix Syop(2); Syop.zero();
  CPPL::zhematrix Szop(2); Szop.zero();

  Sxop(0,1)=0.5;
  Syop(0,1)=-0.5*im;
  Szop(0,0)= 0.5;
  Szop(1,1)=-0.5;

  param pr;
  // if(argc!=3) 
  // {
  //   cerr << "error" << endl;
  //   exit(0);
  // }

  double p = atof(argv[1]);

  pr.K= sin(p);
  pr.J= cos(p);


  // pr.hx = 0.2*0.57735;
  // pr.hy = 0.2*0.57735;
  // pr.hx = 0.12;
  // pr.hy = 0.12;
  double h = atof(argv[2]);
  pr.hx = atof(argv[3]);
  pr.hy = atof(argv[4]);
  pr.hz = atof(argv[5]);




  int N=6;
  vector<double> init(N);
  vector<double> ave(N);

  vector<double> wA,wB;
  vector<CPPL::zcovector> vA;
  vector<CPPL::zcovector> vB;


  double ene_min=5000;

  random_device srd;
  RND rnd(srd());

  int Nr = 100;
  for(int r=0;r<Nr;r++){

    for(int j=0;j<init.size();j++){
      //init[j] = 0.288675;
      init[j] = rnd.d(2.0)-1.0;
    }
    meanfield(ene_min, pr, init, ave, Sxop, Syop, Szop,wA,wB,vA,vB);
  }
  init[0]= 0.499851;
  init[1]= -0.499847;
  init[2]= 0.00545655;
  init[3]= 0.00553761;
  init[4]= -0.0109135;
  init[5]= -0.0110756;
  meanfield(ene_min, pr, init, ave, Sxop, Syop, Szop,wA,wB,vA,vB);

  cout << "Spin # " << ave[0] << " " << ave[2] << " " << ave[4] << " " << ave[1] << " " << ave[3] << " " << ave[5] << endl;


  cout << "MH ene# ";
  cout << boost::format(" %10.3e") % (ene_min/2.0) << endl ;

}
