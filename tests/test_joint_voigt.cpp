// Exact model transforms test the constrained solver independently of sampling.
#include <vector>
#include <iostream>
#include <cstdlib>
#include "../bunch/core/src/joint_voigt_fit.hpp"
struct Context {size_t n=10000;double sigma_p=1.;std::vector<double> t,s;};
struct Moments {std::vector<double> kap,w1,w2,sS;};
int main() {
  Context c;Moments m;constexpr int count=61;
  for(int j=0;j<count;++j){c.t.push_back(.2+.05*j);c.s.push_back(.01);m.kap.push_back(.2);m.w1.push_back(20000.);m.w2.push_back(18000.);m.sS.push_back(.011);}
  for(auto ab:{std::array<double,2>{.3,.6},std::array<double,2>{0.,.7},std::array<double,2>{.4,0.}})
    for(double jump:{-.08,0.,.08}) {
      std::vector<double> C(count),S(count);
      for(int j=0;j<count;++j){double u=c.t[j],phi=std::exp(-ab[0]*u-ab[1]*u*u/2);C[j]=.23*phi;S[j]=jump*phi/u;}
      auto f=sbjointvoigt::fit(c,m,C.data(),S.data(),0,count,nullptr,true);
      if(!f.valid || std::abs(f.B-.23)>1e-7 || std::abs(f.J-jump)>1e-7 || std::abs(f.alpha-ab[0])>1e-6 || std::abs(f.beta-ab[1])>1e-6) {
        std::cerr<<"failed: "<<ab[0]<<" "<<ab[1]<<" "<<jump<<" -> "<<f.valid<<" "<<f.B<<" "<<f.J<<" "<<f.alpha<<" "<<f.beta<<"\n";return 1;
      }
      if(!(f.seB>0 && f.seEq>0)){std::cerr<<"invalid standard error\n";return 1;}
    }
  std::cout<<"9 exact joint Voigt cases (including both boundaries and zero/positive/negative jump) PASS\n";
}
