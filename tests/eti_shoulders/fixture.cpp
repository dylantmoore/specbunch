#include "specbunch.h"
#include "npmle_basis.hpp"
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>
double value(sb_result*r,const char*k){double x;assert(sb_result_scalar(r,k,&x)==0);return x;}
int main(){
 sbnpmle::SplineSide side({0,.1,.4,1},1);std::vector<double> p(side.count),c(side.count);
 for(double x:{0.,.05,.1,.25,.4,.7,1.}){side.latent(x,p.data(),c.data());double sum=0;for(int j=0;j<side.count;++j)sum+=p[j]*side.area[j];assert(std::abs(sum-1)<1e-14);}
 side.latent(1,p.data(),c.data());for(auto x:c)assert(std::abs(x-1)<1e-14);
 // Independent Simpson integration checks the linear Gaussian PDF and CDF.
 for(double x:{-.3,0.,.13,.7,1.3}){
  std::vector<double> pdf(side.count),cdf(side.count),refp(side.count),refc(side.count);
  side.gaussian(x,.07,pdf.data(),cdf.data());
  for(int interval=0;interval<side.intervals;++interval){
   constexpr int steps=4096;double lo=side.knots[interval],h=(side.knots[interval+1]-lo)/steps;
   for(int i=0;i<=steps;++i){double v=lo+i*h,w=(i==0||i==steps)?1:(i%2?4:2);side.latent(v,p.data(),c.data());
    for(int j=0;j<side.count;++j){refp[j]+=h*w/3*p[j]*sbnpmle::normal_pdf((x-v)/.07)/.07;refc[j]+=h*w/3*p[j]*sbnpmle::normal_cdf((x-v)/.07);}
   }
  }
  for(int j=0;j<side.count;++j){assert(std::abs(pdf[j]-refp[j])<1e-9);assert(std::abs(cdf[j]-refc[j])<1e-9);}
 }
 std::vector<double> z(400);for(int i=0;i<400;++i)z[i]=(i-200)/256.;sb_data d{z.size(),z.data(),nullptr};
 std::cout<<"basis,left,right,saez,continuous,lambda,edge,weak,inconclusive\n"<<std::setprecision(17);
 for(int b:{0,1}){
  sb_opts_npmle o;sb_opts_npmle_init(&o);o.B=.08;o.sigma=.03;o.reserved=b;o.core_intervals=4;o.grid_points=25;o.lambda_points=3;o.refinement_points=3;o.threads=1;
  sb_result*n=nullptr;assert(sb_npmle(&d,&o,nullptr,0,&n)==0);
  sb_opts_eti e;sb_opts_eti_init(&e);e.B=o.B;e.net_below=1.34;e.net_above=1;
  double l=value(n,"latent_left"),r=value(n,"latent_right");sb_result*out=nullptr;assert(sb_eti_shoulders(&e,l,r,&out)==0);
  std::cout<<(b?"linear":"cubic")<<','<<l<<','<<r<<','<<value(out,"saez.e")<<','<<value(out,"continuous.e")<<','<<value(n,"lambda")<<','<<value(n,"lambda_at_grid_edge")<<','<<value(n,"shoulders_weak")<<','<<value(n,"shoulders_inconclusive")<<'\n';sb_result_free(out);sb_result_free(n);
  assert(sb_eti_shoulders(&e,0,0,&out)==0);assert(sb_result_status(out)!=0);sb_result_free(out);
  assert(sb_eti_shoulders(&e,-.1,.2,&out)==0);assert(std::isfinite(value(out,"saez.e")));assert(std::isnan(value(out,"continuous.e")));sb_result_free(out);
 }
}
