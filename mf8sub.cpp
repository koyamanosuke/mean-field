#include <iostream>
#include <random>
#include <complex>
#include<fstream>
#include<iomanip>

//#include <boost/format.hpp>

#include <cpplapack/cpplapack.h>
#define Nitr 700000
#define mix 0.01

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
               vector<double>& wE, vector<double>& wF,vector<double>& wG, vector<double>& wH,
               vector<CPPL::zcovector>& vA, vector<CPPL::zcovector>& vB,vector<CPPL::zcovector>& vC, vector<CPPL::zcovector>& vD,
               vector<CPPL::zcovector>& vE, vector<CPPL::zcovector>& vF,vector<CPPL::zcovector>& vG, vector<CPPL::zcovector>& vH)
{

  mf SxA(Sxop,init[0]);
  mf SxB(Sxop,init[1]);
  mf SxC(Sxop,init[2]);
  mf SxD(Sxop,init[3]);
  mf SxE(Sxop,init[4]);
  mf SxF(Sxop,init[5]);
  mf SxG(Sxop,init[6]);
  mf SxH(Sxop,init[7]);
  mf SyA(Syop,init[8]);
  mf SyB(Syop,init[9]);
  mf SyC(Syop,init[10]);
  mf SyD(Syop,init[11]);
  mf SyE(Syop,init[12]);
  mf SyF(Syop,init[13]);
  mf SyG(Syop,init[14]);
  mf SyH(Syop,init[15]);
  mf SzA(Szop,init[16]);
  mf SzB(Szop,init[17]);
  mf SzC(Szop,init[18]);
  mf SzD(Szop,init[19]);
  mf SzE(Szop,init[20]);
  mf SzF(Szop,init[21]);
  mf SzG(Szop,init[22]);
  mf SzH(Szop,init[23]);


  


  bool conv=0;
  long double eps;

  double eneA,eneB,eneC,eneD,eneE,eneF,eneG,eneH;
  for(int i=0;i<Nitr;i++)
    {
      wA.clear();
      wB.clear();
      wC.clear();
      wD.clear();
      wE.clear();
      wF.clear();
      wG.clear();
      wH.clear();
      vA.clear();
      vB.clear();
      vC.clear();
      vD.clear();
      vE.clear();
      vF.clear();
      vG.clear();
      vH.clear();




      CPPL::zhematrix HA = ((2*pr.K+pr.J)*SxB.ave()+pr.J*(SxF.ave()+SxH.ave())-pr.h)*SxA.op()+((2*pr.K+pr.J)*SyF.ave()+pr.J*(SyB.ave()+SyH.ave())-pr.h)*SyA.op()
                          +((2*pr.K+pr.J)*SzH.ave()+pr.J*(SzB.ave()+SzF.ave())-pr.h)*SzA.op();

      CPPL::zhematrix HB = ((2*pr.K+pr.J)*SxA.ave()+pr.J*(SxG.ave()+SxC.ave())-pr.h)*SxB.op()+((2*pr.K+pr.J)*SyG.ave()+pr.J*(SyA.ave()+SyC.ave())-pr.h)*SyB.op()
                          +((2*pr.K+pr.J)*SzC.ave()+pr.J*(SzA.ave()+SzG.ave())-pr.h)*SzB.op();
    
      CPPL::zhematrix HC = ((2*pr.K+pr.J)*SxH.ave()+pr.J*(SxD.ave()+SxB.ave())-pr.h)*SxC.op()+((2*pr.K+pr.J)*SyD.ave()+pr.J*(SyH.ave()+SyB.ave())-pr.h)*SyC.op()
                          +((2*pr.K+pr.J)*SzB.ave()+pr.J*(SzH.ave()+SzD.ave())-pr.h)*SzC.op();

      CPPL::zhematrix HD = ((2*pr.K+pr.J)*SxE.ave()+pr.J*(SxC.ave()+SxG.ave())-pr.h)*SxD.op()+((2*pr.K+pr.J)*SyC.ave()+pr.J*(SyE.ave()+SyG.ave())-pr.h)*SyD.op()
                          +((2*pr.K+pr.J)*SzG.ave()+pr.J*(SzE.ave()+SzC.ave())-pr.h)*SzD.op();

      CPPL::zhematrix HE = ((2*pr.K+pr.J)*SxD.ave()+pr.J*(SxH.ave()+SxF.ave())-pr.h)*SxE.op()+((2*pr.K+pr.J)*SyH.ave()+pr.J*(SyD.ave()+SyF.ave())-pr.h)*SyE.op()
                          +((2*pr.K+pr.J)*SzF.ave()+pr.J*(SzD.ave()+SzH.ave())-pr.h)*SzE.op();

      CPPL::zhematrix HF = ((2*pr.K+pr.J)*SxG.ave()+pr.J*(SxA.ave()+SxE.ave())-pr.h)*SxF.op()+((2*pr.K+pr.J)*SyA.ave()+pr.J*(SyG.ave()+SyE.ave())-pr.h)*SyF.op()
                          +((2*pr.K+pr.J)*SzE.ave()+pr.J*(SzG.ave()+SzA.ave())-pr.h)*SzF.op();

      CPPL::zhematrix HG = ((2*pr.K+pr.J)*SxF.ave()+pr.J*(SxB.ave()+SxD.ave())-pr.h)*SxG.op()+((2*pr.K+pr.J)*SyB.ave()+pr.J*(SyF.ave()+SyD.ave())-pr.h)*SyG.op()
                          +((2*pr.K+pr.J)*SzD.ave()+pr.J*(SzF.ave()+SzB.ave())-pr.h)*SzG.op();

      CPPL::zhematrix HH = ((2*pr.K+pr.J)*SxC.ave()+pr.J*(SxE.ave()+SxA.ave())-pr.h)*SxH.op()+((2*pr.K+pr.J)*SyE.ave()+pr.J*(SyC.ave()+SyA.ave())-pr.h)*SyH.op()
                          +((2*pr.K+pr.J)*SzA.ave()+pr.J*(SzC.ave()+SzE.ave())-pr.h)*SzH.op();


      HA.zheev(wA,vA);
      HB.zheev(wB,vB);
      HC.zheev(wC,vC);
      HD.zheev(wD,vD);
      HE.zheev(wE,vE);
      HF.zheev(wF,vF);
      HG.zheev(wG,vG);
      HH.zheev(wH,vH);
      eneA = wA[0];
      eneB = wB[0];
      eneC = wC[0];
      eneD = wD[0];
      eneE = wE[0];
      eneF = wF[0];
      eneG = wG[0];
      eneH = wH[0];

      eps = 0;
      eps += SxA.calc(vA[0]);
      eps += SxB.calc(vB[0]);
      eps += SxC.calc(vC[0]);
      eps += SxD.calc(vD[0]);
      eps += SxE.calc(vE[0]);
      eps += SxF.calc(vF[0]);
      eps += SxG.calc(vG[0]);
      eps += SxH.calc(vH[0]);
      eps += SyA.calc(vA[0]);
      eps += SyB.calc(vB[0]);
      eps += SyC.calc(vC[0]);
      eps += SyD.calc(vD[0]);
      eps += SyE.calc(vE[0]);
      eps += SyF.calc(vF[0]);
      eps += SyG.calc(vG[0]);
      eps += SyH.calc(vH[0]);
      eps += SzA.calc(vA[0]);
      eps += SzB.calc(vB[0]);
      eps += SzC.calc(vC[0]);
      eps += SzD.calc(vD[0]);
      eps += SzE.calc(vE[0]);
      eps += SzF.calc(vF[0]);
      eps += SzG.calc(vG[0]);
      eps += SzH.calc(vH[0]);



      if(eps < 1e-12) {
          conv = 1;
          break;}
    }

  double total_ene =(eneA + eneB + eneC + eneD +eneE + eneF + eneG + eneH 
                     - 2*pr.K*(SxA.ave()*SxB.ave()+SyA.ave()*SyF.ave()+SzA.ave()*SzH.ave()+SyB.ave()*SyG.ave()+SzB.ave()*SzC.ave()+SxC.ave()*SxH.ave()
                              +SyC.ave()*SyD.ave()+SxD.ave()*SxE.ave()+SzD.ave()*SzG.ave()+SyE.ave()*SyH.ave()+SzE.ave()*SzF.ave()+SxF.ave()*SxG.ave())
                       - pr.J*(SxA.ave()*(SxB.ave()+SxF.ave()+SxH.ave())+SyA.ave()*(SyB.ave()+SyF.ave()+SyH.ave())+SzA.ave()*(SzB.ave()+SzF.ave()+SzH.ave())
                              +SxB.ave()*(SxG.ave()+SxC.ave())+SyB.ave()*(SyG.ave()+SyC.ave())+SzB.ave()*(SzG.ave()+SzC.ave())  
                              +SxC.ave()*(SxH.ave()+SxD.ave())+SyC.ave()*(SyH.ave()+SyD.ave())+SzC.ave()*(SzH.ave()+SzD.ave())   
                              +SxD.ave()*(SxE.ave()+SxG.ave())+SyD.ave()*(SyE.ave()+SyG.ave())+SzD.ave()*(SzE.ave()+SzG.ave())
                              +SxE.ave()*(SxH.ave()+SxF.ave())+SyE.ave()*(SyH.ave()+SyF.ave())+SzE.ave()*(SzH.ave()+SzF.ave())
                              +SxF.ave()*SxG.ave()+SyF.ave()*SyG.ave()+SzF.ave()*SzG.ave()))/8.;     


  // double total_ene =(eneA + eneB + eneC + eneD 
  //                    - 2*pr.K*(SxA.ave()*SxD.ave()+SyA.ave()*SyD.ave()+SzA.ave()*SzB.ave()+SxB.ave()*SxC.ave()+SyB.ave()*SyC.ave()+SzC.ave()*SzD.ave())
  //                    - 2*pr.J*(SxA.ave()*SxD.ave()+SyA.ave()*SyD.ave()+SzA.ave()*SzD.ave()+SxB.ave()*SxC.ave()+SyB.ave()*SyC.ave()+SzB.ave()*SzC.ave())
  //                    -   pr.J*(SxA.ave()*SxB.ave()+SyA.ave()*SyB.ave()+SzA.ave()*SzB.ave()+SxC.ave()*SxD.ave()+SyC.ave()*SyD.ave()+SzC.ave()*SzD.ave()))/2.;            


  if(conv){
    if(total_ene < ene_min){
      ene_min = total_ene;
      ave[0] =  SxA.ave();
      ave[1] =  SxB.ave();
      ave[2] =  SxC.ave();
      ave[3] =  SxD.ave();
      ave[4] =  SxE.ave();
      ave[5] =  SxF.ave();
      ave[6] =  SxG.ave();
      ave[7] =  SxH.ave();
      ave[8] =  SyA.ave();
      ave[9] =  SyB.ave();
      ave[10] = SyC.ave();
      ave[11] = SyD.ave();
      ave[12] = SyE.ave();
      ave[13] = SyF.ave();
      ave[14] = SyG.ave();
      ave[15] = SyH.ave();
      ave[16] = SzA.ave();
      ave[17] = SzB.ave();
      ave[18] = SzC.ave();
      ave[19] = SzD.ave();
      ave[20] = SzE.ave();
      ave[21] = SzF.ave();
      ave[22] = SzG.ave();
      ave[23] = SzH.ave();

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

  ofstream ofs("KH8sub111.txt");
  ofstream ofs1("KH8data111.txt");

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

  int N=24;
  vector<long double> init(N);
  vector<long double> ave(N);

  vector<double> wA,wB,wC,wD,wE,wF,wG,wH;
  vector<CPPL::zcovector> vA;
  vector<CPPL::zcovector> vB;
  vector<CPPL::zcovector> vC;
  vector<CPPL::zcovector> vD;
  vector<CPPL::zcovector> vE;
  vector<CPPL::zcovector> vF;
  vector<CPPL::zcovector> vG;
  vector<CPPL::zcovector> vH;

  

  //for(double p=-3.14153;p<3.1415;p+=0.05){
  for(double p=-1.56 ;p<-0.489;p+=0.05){
    //double p = 0.75;
    //double p=atof(argv[1]);
    double ene_min=5000;
    param pr;
    //pr.h=atof(argv[2]);
    pr.K = sin(p);
    pr.J = cos(p);
    //pr.h =0;
    cout << "                \x1b[31m θ = " << p << "\x1b[m" << endl;


    //random_device srd;
    //RND rnd(srd());

    for(pr.h=0;pr.h<6;pr.h+=0.005){
      int Nr =  1000;
      // for(int i=0;i<Nr;i++){
        
      //   for(int j=0;j<init.size();j++){
            
      //     init[j]= rnd.d(2.0)-1.0;

      //   }
      //   meanfield(ene_min, pr, init, ave, Sxop, Syop, Szop,wA,wB,wC,wD,wE,wF,wG,wH,
      //             vA,vB,vC,vD,vE,vF,vG,vH);
      // }
    
    if(p<-0.489){
      init[0] = -0.22965;
      init[1] = -0.229652;
      init[2] = 0.288675;    
      init[3] = -0.229652;
      init[4] = -0.229652;
      init[5] = 0.288675;
      init[6] = 0.380158;
      init[7] = 0.380158;
      init[8] = 0.380158;
      init[9] = -0.229652;
      init[10]= 0.288675;
      init[11]= 0.380158;
      init[12]= -0.229652;
      init[13]= 0.288675;
      init[14]= -0.229652;
      init[15]= -0.229652;
      init[16]= -0.229652;
      init[17]= 0.380158;
      init[18]= 0.288675;
      init[19]= -0.229652;
      init[20]= 0.380158;
      init[21]= 0.288675;
      init[22]= -0.229652;
      init[23]= -0.229652;
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

    meanfield(ene_min, pr, init, ave, Sxop, Syop, Szop,wA,wB,wC,wD,wE,wF,wG,wH,
                  vA,vB,vC,vD,vE,vF,vG,vH);


      cout << "# ";
      cout << "h:" << pr.h << "  ene_min:" << setprecision(4) << ene_min << endl;
      for(int i=0;i<N;i++){
        cout << setprecision(7) <<  ave[i] << "  ";
      }
      
      long double diff = 0.;
      diff += fabsl(ave[0]-0.28867);
      diff += fabsl(ave[1]-0.28867);
      diff += fabsl(ave[2]-0.28867);
      diff += fabsl(ave[3]-0.28867);
      diff += fabsl(ave[4]-0.28867);
      diff += fabsl(ave[5]-0.28867);
      diff += fabsl(ave[6]-0.28867);
      diff += fabsl(ave[7]-0.28867);
      diff += fabsl(ave[8]-0.28867);
      diff += fabsl(ave[9]-0.28867);
      diff += fabsl(ave[10]-0.28867);
      diff += fabsl(ave[11]-0.28867);
      diff += fabsl(ave[12]-0.28867);
      diff += fabsl(ave[13]-0.28867);
      diff += fabsl(ave[14]-0.28867);
      diff += fabsl(ave[15]-0.28867);
      diff += fabsl(ave[16]-0.28867);
      diff += fabsl(ave[17]-0.28867);
      diff += fabsl(ave[18]-0.28867);
      diff += fabsl(ave[19]-0.28867);
      diff += fabsl(ave[20]-0.28867);
      diff += fabsl(ave[21]-0.28867);
      diff += fabsl(ave[22]-0.28867);
      diff += fabsl(ave[23]-0.28867);
       if(diff<1e-3){
         break;
       }
      cout << endl << endl;
      ofs1 << p << " " << pr.h << " " << ave[0] << " " << ave[1] << " " << ave[2] << " " << ave[3] << " " << ave[4] << " " << ave[5] << " " 
           << ave[6] << " " << ave[7] << " " << ave[8] << " " << ave[9] << " " << ave[10] << " " << ave[11]
           << ave[12] << " " << ave[13] << " " << ave[14] << " " << ave[15] << " " << ave[16] << " " << ave[17] << " " 
           << ave[18] << " " << ave[19] << " " << ave[20] << " " << ave[21] << " " << ave[22] << " " << ave[23] << endl;
    }
    
    ofs << p << " " << pr.h << endl;
    cout << endl << "--------------------------------------------------------------------------------------------------------------------------------" << endl;
    ave.clear();
  }
  
}