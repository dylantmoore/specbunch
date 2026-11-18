#pragma once
#include <algorithm>
#include <cmath>
#include <vector>
#include "sbc_internal.hpp"
#include "specbunch.h"
namespace sbtransform {
// Diagnostic evaluation only: never feeds the estimator or window selection.
// Fixed observation chunks and fixed-order reductions make thread count irrelevant.
inline void add(const sb_data& data, int threads, sbc::Result& r) {
  const bool joint=r.strings["method"]=="joint";
  const char* columns[]={"t","t_scaled","cosine","cosine_se","cosine_fit","sine_scaled","sine_scaled_se","sine_fit"};
  for(int j=0;j<(joint?8:5);++j)r.setv(std::string("transform.")+columns[j],{});
  r.set("n_points",NAN);
  auto &chosen=r.vectors["windows.selected"];
  for(size_t i=0;i<chosen.size();++i)if(chosen[i]!=0)r.set("n_points",r.vectors["windows.n_points"][i]);
  const double pilot=r.scalars["sigma_pilot"];
  const auto& lo=r.vectors["windows.window_lo"];
  const auto& hi=r.vectors["windows.window_hi"];
  if(!(pilot>0 && std::isfinite(pilot)) || lo.empty() || hi.empty())return;
  const double first=std::max(0.,*std::min_element(lo.begin(),lo.end())-.1/pilot);
  const double last=*std::max_element(hi.begin(),hi.end())+.2/pilot;
  if(!(last>first && std::isfinite(last)))return;
  constexpr int points=401,chunks=64,block=256,cols=6;
  const double step=(last-first)/(points-1);
  std::vector<double> z,w;z.reserve(data.n);w.reserve(data.n);
  double W=0,W2=0;
  for(size_t i=0;i<data.n;++i) {
    const double wi=data.w?data.w[i]:1.;
    if(!std::isfinite(data.z[i]) || !std::isfinite(wi) || wi<=0)continue;
    z.push_back(data.z[i]);w.push_back(wi);W+=wi;W2+=wi*wi;
  }
  if(!(W>0))return;
  std::vector<double> partial((size_t)chunks*points*cols,0.);
  const size_t nblocks=(z.size()+block-1)/block,per=(nblocks+chunks-1)/chunks;
  sbc::ThreadGuard guard(threads);
#pragma omp parallel for schedule(static)
  for(int ch=0;ch<chunks;++ch) {
    double *p=&partial[(size_t)ch*points*cols];
    double cc[block],ss[block],cd[block],sd[block];
    for(size_t b=ch*per;b<std::min(nblocks,(ch+1)*per);++b) {
      const size_t start=b*block;const int count=(int)std::min((size_t)block,z.size()-start);
      for(int i=0;i<count;++i){cd[i]=std::cos(step*z[start+i]);sd[i]=std::sin(step*z[start+i]);}
      for(int j=0;j<points;++j) {
        if(j%32==0)for(int i=0;i<count;++i){const double a=(first+step*j)*z[start+i];cc[i]=std::cos(a);ss[i]=std::sin(a);}
        else for(int i=0;i<count;++i){double c=cc[i];cc[i]=c*cd[i]-ss[i]*sd[i];ss[i]=ss[i]*cd[i]+c*sd[i];}
        double sums[cols]={};
        for(int i=0;i<count;++i) {
          double wi=w[start+i],sq=wi*wi,c=cc[i],s=ss[i];
          sums[0]+=wi*c;sums[1]+=sq*c;sums[2]+=sq*c*c;
          if(joint){sums[3]+=wi*s;sums[4]+=sq*s;sums[5]+=sq*s*s;}
        }
        for(int k=0;k<cols;++k)p[j*cols+k]+=sums[k];
      }
    }
  }
  std::vector<double> values[8];for(auto &v:values)v.resize(points);
  const double B=r.scalars["B"],sigma=r.scalars["sigma"],a=r.strings["friction"]=="voigt"?r.scalars["a"]:0.;
  const double jump=joint?r.scalars["jump"]:0.;
  for(int j=0;j<points;++j) {
    double sums[cols]={};for(int ch=0;ch<chunks;++ch)for(int k=0;k<cols;++k)sums[k]+=partial[((size_t)ch*points+j)*cols+k];
    const double t=first+step*j,c=sums[0]/W,s=sums[3]/W;
    const double phi=std::exp(-a*std::abs(t)-.5*sigma*sigma*t*t);
    values[0][j]=t;values[1][j]=t*pilot;values[2][j]=c;
    values[3][j]=std::sqrt(std::max(0.,sums[2]-2*c*sums[1]+c*c*W2))/W;
    values[4][j]=r.status==0?B*phi:NAN;
    values[5][j]=t*s;values[6][j]=t*std::sqrt(std::max(0.,sums[5]-2*s*sums[4]+s*s*W2))/W;
    values[7][j]=r.status==0?jump*phi:NAN;
  }
  for(int j=0;j<(joint?8:5);++j)r.setv(std::string("transform.")+columns[j],std::move(values[j]));
}
}
