#include <iostream>
#include <random>
#include <complex>
#include<fstream>
#include<iomanip>

//#include <boost/format.hpp>

#include <cpplapack/cpplapack.h>
#define Nitr 300000
#define mix 0.05

using namespace std;

class mf{
private:
  long double val;
  const CPPL::zhematrix& mat;
  long double init;
public:
  mf(const CPPL::zhematrix& mat_): mat(mat_), init(1000){val=init;}
  mf(const CPPL::zhematrix& mat_, const double init_): mat(mat_), init(init_){val=init;}
  CPPL::zhematrix op(){
    return mat;
  }
  double calc(const CPPL::zcovector& v){
    CPPL::zcovector vt = mat*v;
    complex<long double> sum(0,0);
    for(int i=0;i<v.l;i++)
      sum += conj(v(i))*vt(i);

    long double val_old = val;
    val = std::real(sum)*mix + val_old*(1-mix);
    return fabsl(val-val_old);
  }
  double ave(){return val;};
};

struct param{
  double J;
  double K;
  double h;
};
void meanfield(double& ene_min, const param& pr,
               vector<long double>& init, vector<long double>& ave,
               const CPPL::zhematrix& Sxop, const CPPL::zhematrix& Syop, const CPPL::zhematrix& Szop,
               vector<double>& wA, vector<double>& wB,vector<double>& wC, vector<double>& wD,
               vector<CPPL::zcovector>& vA, vector<CPPL::zcovector>& vB,vector<CPPL::zcovector>& vC, vector<CPPL::zcovector>& vD)
{

  mf SxA(Sxop,init[0]);
  mf SxB(Sxop,init[1]);
  mf SxC(Sxop,init[2]);
  mf SxD(Sxop,init[3]);
  mf SyA(Syop,init[4]);
  mf SyB(Syop,init[5]);
  mf SyC(Syop,init[6]);
  mf SyD(Syop,init[7]);
  mf SzA(Szop,init[8]);
  mf SzB(Szop,init[9]);
  mf SzC(Szop,init[10]);
  mf SzD(Szop,init[11]);

  


  bool conv=0;
  long double eps;

  double eneA, eneB,eneC,eneD;
  for(int i=0;i<Nitr;i++)
    {
      wA.clear();
      wB.clear();
      wC.clear();
      wD.clear();
      vA.clear();
      vB.clear();
      vC.clear();
      vD.clear();




      CPPL::zhematrix HA = ((2*pr.K+2*pr.J)*SxD.ave()+pr.J*SxB.ave()-pr.h)*SxA.op()+((2*pr.K+2*pr.J)*SyD.ave()+pr.J*SyB.ave()-pr.h)*SyA.op()
                           +(2*pr.J*SzD.ave()+(2*pr.K+pr.J)*SzB.ave()-pr.h)*SzA.op();

      CPPL::zhematrix HB = ((2*pr.K+2*pr.J)*SxC.ave()+pr.J*SxA.ave()-pr.h)*SxB.op()+((2*pr.K+2*pr.J)*SyC.ave()+pr.J*SyA.ave()-pr.h)*SyB.op()
                           +(2*pr.J*SzC.ave()+(2*pr.K+pr.J)*SzA.ave()-pr.h)*SzB.op();

      CPPL::zhematrix HC = ((2*pr.K+2*pr.J)*SxB.ave()+pr.J*SxD.ave()-pr.h)*SxC.op()+((2*pr.K+2*pr.J)*SyB.ave()+pr.J*SyD.ave()-pr.h)*SyC.op()
                           +(2*pr.J*SzB.ave()+(2*pr.K+pr.J)*SzD.ave()-pr.h)*SzC.op();

      CPPL::zhematrix HD = ((2*pr.K+2*pr.J)*SxA.ave()+pr.J*SxC.ave()-pr.h)*SxD.op()+((2*pr.K+2*pr.J)*SyA.ave()+pr.J*SyC.ave()-pr.h)*SyD.op()
                           +(2*pr.J*SzA.ave()+(2*pr.K+pr.J)*SzC.ave()-pr.h)*SzD.op();
    
                         


//cout << HA << HB << endl;

      HA.zheev(wA,vA);
      HB.zheev(wB,vB);
      HC.zheev(wC,vC);
      HD.zheev(wD,vD);
      eneA = wA[0];
      eneB = wB[0];
      eneC = wC[0];
      eneD = wD[0];
      //cout << eneA << " " << eneB << endl;

      eps = 0;
      eps += SxA.calc(vA[0]);
      eps += SxB.calc(vB[0]);
      eps += SxC.calc(vC[0]);
      eps += SxD.calc(vD[0]);
      eps += SyA.calc(vA[0]);
      eps += SyB.calc(vB[0]);
      eps += SyC.calc(vC[0]);
      eps += SyD.calc(vD[0]);
      eps += SzA.calc(vA[0]);
      eps += SzB.calc(vB[0]);
      eps += SzC.calc(vC[0]);
      eps += SzD.calc(vD[0]);



      if(eps < 1e-10) {
          conv = 1;
          break;}
    }

  double total_ene =(eneA + eneB + eneC + eneD 
                     - 2*pr.K*(SxA.ave()*SxD.ave()+SyA.ave()*SyD.ave()+SzA.ave()*SzB.ave()+SxB.ave()*SxC.ave()+SyB.ave()*SyC.ave()+SzC.ave()*SzD.ave())
                     - 2*pr.J*(SxA.ave()*SxD.ave()+SyA.ave()*SyD.ave()+SzA.ave()*SzD.ave()+SxB.ave()*SxC.ave()+SyB.ave()*SyC.ave()+SzB.ave()*SzC.ave())
                     -   pr.J*(SxA.ave()*SxB.ave()+SyA.ave()*SyB.ave()+SzA.ave()*SzB.ave()+SxC.ave()*SxD.ave()+SyC.ave()*SyD.ave()+SzC.ave()*SzD.ave()))/2.;            

  if(conv){
    if(total_ene < ene_min){
      ene_min = total_ene;
      ave[0] = SxA.ave();
      ave[1] = SxB.ave();
      ave[2] = SxC.ave();
      ave[3] = SxD.ave();
      ave[4] = SyA.ave();
      ave[5] = SyB.ave();
      ave[6] = SyC.ave();
      ave[7] = SyD.ave();
      ave[8] = SzA.ave();
      ave[9] = SzB.ave();
      ave[10] =SzC.ave();
      ave[11] =SzD.ave();

    }
  }

  // cerr << "conv = " << conv << endl;
  // cerr << "eps  = " << eps << endl;
  // cerr << "total energy  = " << total_ene << endl << endl;
}

// void calc_element(const vector<CPPL::zcovector>& vec, const CPPL::zhematrix& mat, vector<complex<double> >& op)
// {
//   op.resize(vec.size()-1);
//   {
//     CPPL::zcovector tv = mat*vec[0];
//     for(int j=1;j<vec.size();j++){
//       complex<double> sum(0,0);
//       for(int i=0;i<vec[j].l;i++){
//         sum+= conj(vec[j](i))*tv(i);
//       }
//       op[j-1]=sum;
//     }
//   }
// }

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

int main()
{

  // ofstream ofs("KH4sub111.txt");
  // ofstream ofs1("KHdata111.txt");

  complex<double> im(0,1);

  CPPL::zhematrix Sxop(2); Sxop.zero();
  CPPL::zhematrix Syop(2); Syop.zero();
  CPPL::zhematrix Szop(2); Szop.zero();

  Sxop(0,1)=0.5;
  Syop(0,1)=-0.5*im;
  Szop(0,0)= 0.5;
  Szop(1,1)=-0.5;

//   param pr;
//   if(argc!=3) exit(0);

//   pr.J=atof(argv[1]);
//   pr.K=atof(argv[2]);
//   double p=atof(argv[1]);
//   pr.h=atof(argv[2]);

//   pr.hx = atof(argv[4]);
//   pr.hy = atof(argv[5]);
//   pr.hz = atof(argv[6]);

  int N=12;
  vector<long double> init(N);
  vector<long double> ave(N);

  vector<double> wA,wB,wC,wD;
  vector<CPPL::zcovector> vA;
  vector<CPPL::zcovector> vB;
  vector<CPPL::zcovector> vC;
  vector<CPPL::zcovector> vD;

  

  //for(double p=-3.14153;p<3.1415;p+=0.05){
  for(double p=-1.5 ;p<1.6;p+=0.01){
    //double p = 0.75;
    double ene_min=5000;
    param pr;
    pr.K = sin(p);
    //pr.J = cos(p);
    pr.J=0.;
    //pr.h =0;
    cout << "                \x1b[31m θ = " << p << "\x1b[m" << endl;


    random_device srd;
    RND rnd(srd());

     //for(pr.h=0;pr.h<8;pr.h+=0.01){
       pr.h=0.;
      // int Nr = 100;
      // for(int i=0;i<Nr;i++){
        
      //   for(int j=0;j<init.size();j++){
            
      //     init[j]= rnd.d(2.0)-1.0;

      //   }
      //   meanfield(ene_min, pr, init, ave, Sxop, Syop, Szop,wA,wB,wC,wD,vA,vB,vC,vD);
      // }
    
    if(p<-0.489){
      init[0] = 2.12225e-08;
      init[1] = -2.12225e-08;
      init[2] = 2.12225e-08;    
      init[3] = -2.12225e-08;
      init[4] = 3.88933e-08;
      init[5] = -3.88933e-08;
      init[6] = 3.88933e-08;
      init[7] = -3.88933e-08;
      init[8] = -0.5;
      init[9] = -0.5;
      init[10]= 0.5;
      init[11]= 0.5;
    }

    else if(p>-0.489 && p<1.570){
      init[0] = 0.0892472;
      init[1] = -0.0892472;
      init[2] = 0.0892472;    
      init[3] = -0.0892472;
      init[4] = 0.130571;
      init[5] = -0.130571;
      init[6] = 0.130571;
      init[7] = -0.130571;
      init[8] = 0.474327;
      init[9] = -0.474327;
      init[10]= 0.474327;
      init[11]= -0.474327;
    }

     else if(p>1.570){
      init[0] = 7.11816e-08;
      init[1] = 7.11817e-08;
      init[2] = -7.11817e-08;    
      init[3] = -7.11817e-08;
      init[4] = -7.0448e-08;
      init[5] = -7.0448e-08;
      init[6] = 7.0448e-08;
      init[7] = 7.0448e-08;
      init[8] = 0.5;
      init[9] = -0.5;
      init[10]= -0.5;
      init[11]= 0.5;
    }

     meanfield(ene_min, pr, init, ave, Sxop, Syop, Szop,wA,wB,wC,wD,vA,vB,vC,vD);


      cout << "# ";
      cout << "h:" << pr.h << "  ene_min:" << setprecision(4) << ene_min << endl;

    //   for(int i=0;i<N;i++){
    //     cout << setprecision(9) <<  ave[i] << "  ";
    //   }
      
    //   long double diff = 0.;
    //   diff += fabsl(ave[0]-0.288675134);
    //   diff += fabsl(ave[1]-0.288675134);
    //   diff += fabsl(ave[2]-0.288675134);
    //   diff += fabsl(ave[3]-0.288675134);
    //   diff += fabsl(ave[4]-0.288675134);
    //   diff += fabsl(ave[5]-0.288675134);
    //   diff += fabsl(ave[6]-0.288675134);
    //   diff += fabsl(ave[7]-0.288675134);
    //   diff += fabsl(ave[8]-0.288675134);
    //   diff += fabsl(ave[9]-0.288675134);
    //   diff += fabsl(ave[10]-0.288675134);
    //   diff += fabsl(ave[11]-0.288675134);
    //    if(diff<5e-7){
    //      break;
    //    }
    //   cout << endl << endl;
    //   // ofs1 << p << " " << pr.h << " " << ave[0] << " " << ave[1] << " " << ave[2] << " " << ave[3] << " " << ave[4] << " " << ave[5] << " " 
    //   //      << ave[6] << " " << ave[7] << " " << ave[8] << " " << ave[9] << " " << ave[10] << " " << ave[11] << endl;
    // //}
    
    // //ofs << p << " " << pr.h << endl;
    // cout << endl << "--------------------------------------------------------------------------------------" << endl;
    // ave.clear();
  }
  
}