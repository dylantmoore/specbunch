#include <algorithm>
#include <cmath>
#include <limits>
// Bracket inversion of the small Gaussian tail: no cancellation near one.
static double invnormal(double p){
 if(p<=0)return -INFINITY;if(p>=1)return INFINITY;
 double tail=p<.5?p:1-p,lo=0,hi=40;
 for(int j=0;j<64;++j){double mid=(lo+hi)/2;if(.5*std::erfc(mid/std::sqrt(2.))>tail)lo=mid;else hi=mid;}
 return (p<.5?-1:1)*(lo+hi)/2;
}
extern "C" double normal_quantile_test(double p){return invnormal(p);}
extern "C" int probit_screen(const double*x,const double*f,int n,double B,double sigma,double*margin,int*endpoint){
 if(n<=0||!std::isfinite(B)||B<0||B>1||!std::isfinite(sigma)||sigma<0)return -1;
 double mass=1-B,eps=std::sqrt(std::log(2/.001)/(2*n)),rl=0,ru=1,best=-1;int hit=-1;
 constexpr double pad=2e-12;
 for(int j=0;j<66;++j){
  double G=sigma>0?.5*std::erfc(-x[j]/(sigma*std::sqrt(2.))):(x[j]>=0?1.:0.);
  double L=std::max(0.,f[j]-eps-pad),U=std::min(1.,f[j]+eps+pad);
  double point=std::max(B*G-U,L-(B*G+mass));if(point>best){best=point;hit=j;}
  if(mass==0)continue;
  double l=std::max(0.,(L-B*G)/mass),u=std::min(1.,(U-B*G)/mass);
  if(l>u){double gap=mass*(l-u);if(gap>best){best=gap;hit=j;}continue;}
  l=std::max(0.,l-pad);u=std::min(1.,u+pad);
  if(j==0){rl=l;ru=u;}
  else{
   if(sigma>0){double delta=(x[j]-x[j-1])/sigma;
    double reachable=ru<=0?0:(ru>=1?1:.5*std::erfc(-(invnormal(ru)+delta)/std::sqrt(2.)));
    ru=std::min(u,std::min(1.,reachable+pad));
   }else ru=u; // sigma=0 imposes monotonicity only; do not use g_delta(0)=0.
   rl=std::max(rl,l);
  }
  double gap=mass*(rl-ru);if(gap>best){best=gap;hit=j;}
 }
 *margin=best;*endpoint=hit;return best>1e-10;
}
