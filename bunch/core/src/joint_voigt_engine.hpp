#pragma once
#include "joint_voigt_fit.hpp"
namespace sbjointvoigt {
struct Table {int M=0;std::vector<Fit> fits;std::vector<double> p;std::vector<int> rank;};
inline Table windows(const sbjoint_reference::Ctx& c,const sbjoint_reference::Mom& mom,const std::vector<double>& Yc,const std::vector<double>& Ys,const std::vector<double>& signs) {
  using std::vector;const int ngrid=c.M;Table table;table.M=ngrid;
  table.fits.resize((size_t)ngrid*ngrid);table.p.assign(table.fits.size(),NAN);
  vector<double> cosine_model((size_t)ngrid*ngrid*ngrid,0.),sine_model(cosine_model.size(),0.);
  vector<int> all;
  for(int first=0;first<ngrid;++first)for(int last=first+sbjoint_reference::MINPTS-1;last<ngrid;++last)all.push_back(first*ngrid+last);
#pragma omp parallel for schedule(dynamic,4)
  for(size_t k=0;k<all.size();++k) {
    int id=all[k],first=id/ngrid,last=id%ngrid;
    Fit f=fit(c,mom,c.C.data(),mom.S.data(),first,last-first+1,nullptr,true);
    table.fits[id]=f;
    if(!f.valid)continue;
    for(int j=first;j<=last;++j) {
      double u=c.t[j]*c.sigma_p,e=std::exp(-f.alpha*u-.5*f.beta*u*u);
      cosine_model[(size_t)id*ngrid+j]=f.B*e;sine_model[(size_t)id*ngrid+j]=f.J*e/u;
    }
  }
  vector<int> ids;for(int id:all)if(table.fits[id].valid&&std::isfinite(table.fits[id].D))ids.push_back(id);
  vector<int> exceed(ids.size(),0);
#pragma omp parallel
  {
    vector<int> local(ids.size(),0);double C[sbjoint_reference::NGPAD],S[sbjoint_reference::NGPAD],dc[sbjoint_reference::NGPAD],ds[sbjoint_reference::NGPAD];
#pragma omp for schedule(dynamic,4)
    for(int b=0;b<sbjoint_reference::BOOT_B;++b) {
      for(int j=0;j<ngrid;++j){dc[j]=(Yc[(size_t)b*sbjoint_reference::NGPAD+j]-c.C[j]*signs[b])/c.n;ds[j]=(Ys[(size_t)b*sbjoint_reference::NGPAD+j]-mom.S[j]*signs[b])/c.n;}
      for(size_t k=0;k<ids.size();++k) {
        int id=ids[k],first=id/ngrid,last=id%ngrid;
        for(int j=first;j<=last;++j){C[j]=cosine_model[(size_t)id*ngrid+j]+dc[j];S[j]=sine_model[(size_t)id*ngrid+j]+ds[j];}
        Fit f=fit(c,mom,C,S,first,last-first+1,&table.fits[id],false);
        if(f.valid && f.D>=table.fits[id].D)++local[k];
      }
    }
#pragma omp critical
    for(size_t k=0;k<ids.size();++k)exceed[k]+=local[k]; // integer reduction is exact
  }
  for(size_t k=0;k<ids.size();++k)table.p[ids[k]]=(1.+exceed[k])/(sbjoint_reference::BOOT_B+1.);
  table.rank=ids;
  std::stable_sort(table.rank.begin(),table.rank.end(),[&](int x,int y) {
    bool passx=table.p[x]>=.10,passy=table.p[y]>=.10;if(passx!=passy)return passx;
    if(!passx && table.p[x]!=table.p[y])return table.p[x]>table.p[y];
    int nx=x%ngrid-x/ngrid,ny=y%ngrid-y/ngrid;if(nx!=ny)return nx>ny;
    return x/ngrid<y/ngrid;
  });
  return table;
}
struct Output {double pilot=NAN;std::string reason;std::vector<double> table,parameters,u;};
inline Output run(const std::vector<double>& z,uint64_t seed) {
  using namespace sbjoint_reference;Output out;Pilot pilot=ensemble_sigma(z);out.pilot=pilot.sigma_p;
  if(!(std::isfinite(out.pilot)&&out.pilot>0)){out.reason="no_frequency_band";return out;}
  double start=UGRID_LO/out.pilot,delta=UGRID_DU/out.pilot;
  vector<double> C,cosine_second_moment,S,sine_second_moment,CS;cf_recurrence(z,start,delta,NGMAX,C,cosine_second_moment);sin_moments(z,start,delta,NGMAX,S,sine_second_moment,CS);
  Ctx initial=make_ctx(z.size(),out.pilot,NG0,C,cosine_second_moment);Mom initial_mom=make_mom(initial,S,sine_second_moment,CS);
  vector<double> Yc((size_t)BOOT_B*NGPAD,0.),Ys(Yc.size(),0.),signs,temp;
  boot_g0<false>(z,start,delta,seed,temp,&signs);
  for(int b=0;b<BOOT_B;++b)std::copy_n(&temp[(size_t)b*MP0],NG0,&Yc[(size_t)b*NGPAD]);
  boot_g0<true>(z,start,delta,seed,temp,nullptr);
  for(int b=0;b<BOOT_B;++b)std::copy_n(&temp[(size_t)b*MP0],NG0,&Ys[(size_t)b*NGPAD]);
  sbjointvoigt::Table ranked=windows(initial,initial_mom,Yc,Ys,signs);
  int size=NG0;
  if(!ranked.rank.empty()) {
    const Fit& f=ranked.fits[ranked.rank[0]];
    for(int j=NG0-1;j<NGMAX;++j) {
      double u=UGRID_LO+j*UGRID_DU,t=u/out.pilot,e=std::exp(-f.alpha*u-.5*f.beta*u*u);
      double mc=f.B*e,ms=f.J*e/u;
      double noise=std::sqrt(std::max(0.,cosine_second_moment[j]-C[j]*C[j]+sine_second_moment[j]-S[j]*S[j])/z.size());
      if(std::hypot(mc,ms)>=noise)size=j+1;
      (void)t;
    }
  }
  if(size>NG0) {
    int width=std::min(NGPAD,((std::max(NG1,size)+15)/16)*16);
    boot_cols<false>(z,start,delta,seed,64,width,temp,nullptr);
    for(int b=0;b<BOOT_B;++b)std::copy(&temp[(size_t)b*NGPAD+NG0],&temp[(size_t)b*NGPAD+width],&Yc[(size_t)b*NGPAD+NG0]);
    temp.clear();boot_cols<true>(z,start,delta,seed,64,width,temp,nullptr);
    for(int b=0;b<BOOT_B;++b)std::copy(&temp[(size_t)b*NGPAD+NG0],&temp[(size_t)b*NGPAD+width],&Ys[(size_t)b*NGPAD+NG0]);
    Ctx c=make_ctx(z.size(),out.pilot,size,C,cosine_second_moment);Mom mom=make_mom(c,S,sine_second_moment,CS);ranked=windows(c,mom,Yc,Ys,signs);out.u=c.u;
  }else out.u=initial.u;
  if(ranked.rank.empty())out.reason="no_coherent_window";
  for(int id:ranked.rank) {
    const Fit& f=ranked.fits[id];double a=f.alpha*out.pilot,s2=f.beta*out.pilot*out.pilot,eq=a+std::sqrt(a*a+s2);
    double row[10]={(double)(id/size),(double)(id%size),f.B,eq,f.seB,f.seEq,ranked.p[id],f.valid?1.:0.,NAN,f.J/out.pilot};
    out.table.insert(out.table.end(),row,row+10);out.parameters.push_back(a);out.parameters.push_back(s2);
  }
  return out;
}
}
