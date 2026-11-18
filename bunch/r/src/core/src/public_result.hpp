#pragma once
#include <cmath>
#include <string>
#include <map>
#include "sbc_internal.hpp"
#include "specbunch.h"
namespace sbpublic {
inline void finish(sbc::Result &r,const std::string &method) {
  const double nan=std::numeric_limits<double>::quiet_NaN();
  auto get=[&](const std::string &k){auto i=r.scalars.find(k);return i==r.scalars.end()?nan:i->second;};
  const bool mean=method=="buncher_mean";
  auto old=r.scalars;auto vec=r.vectors;
  std::map<std::string,double> out;
  for(const std::string k:{"n","n_input","n_eff","weighted","W","seed","threads_used","B","sigma","sigma_eq","a_voigt","s2_voigt","tier","rank","flags","mu","se","sigma_pilot","band_lo","band_hi"})
    if(old.count(k))out[k]=old[k];
  for(auto kv:std::map<std::string,std::string>{{"se_logB","se_log_B"},{"se_logS","se_log_sigma"},{"se_log_sigma_eq","se_log_sigma_eq"},{"p","fit_p"},{"aJ","jump"},{"pw_median","pointwise_median"},{"sigma_hat","sigma_pilot"}})
    if(old.count(kv.first))out[kv.second]=old[kv.first];
  if(!mean)for(auto type:{"ols","lad"}) {
    std::string prefix=std::string(type)=="ols"?"loglinear_ols.":"loglinear_lad.";
    for(auto k:{"B","sigma","status","m","ndrop"}) if(old.count(prefix+k))out[std::string("loglinear_")+type+"."+k]=old[prefix+k];
  }
  const char* flags[]={"cdf_inconsistent","moments_inconsistent","no_valid_distribution","distribution_check_inconclusive","residual_cdf_inconsistent","poor_fit"};
  int bits=std::isfinite(get("flags"))?(int)get("flags"):0;std::string flagtext;
  for(int j=0;j<6;++j) {out[flags[j]]=(bits&(1<<j))?1:0;if(bits&(1<<j))flagtext+=(flagtext.empty()?"":",")+std::string(flags[j]);}
  out["fit_converged"]=r.status==0;out["coherent_values"]=r.status==0;
  std::string reason;
  if(r.status!=0) {
    auto internal=r.strings["fit_failure"];
    reason=internal=="no_frequency_band" || internal.find("pilot")!=std::string::npos || (std::isfinite(get("sigma_pilot"))&&get("sigma_pilot")<=0)?"no_frequency_band":"no_coherent_window";
    if(r.strings["reason"]=="insufficient_data")reason="insufficient_data";
    r.status=reason=="insufficient_data"?1:reason=="no_frequency_band"?2:3;
  }
  out["status"]=r.status;
  std::string policy=get("policy_mode")==2?"first":"tiered";
  r.strings={{"method",mean?method:((method=="gaussian_jump" || method=="joint_voigt")?"joint":"cosine")},{"friction",mean?"":((method=="voigt" || method=="joint_voigt")?"voigt":"gaussian")},{"estimator",mean?"cosine-ratio fit":((method=="gaussian_jump" || method=="joint_voigt")?"Weighted nonlinear least squares of the empirical cosine and sine transforms in levels, with fixed inverse-covariance weights from their empirical covariance, solved by Gauss-Newton/Newton.":"Weighted nonlinear least squares of the empirical cosine transform in levels, with fixed inverse-variance weights from the transform variance, solved by Gauss-Newton/Newton.")},{"policy",policy},{"reason",reason},{"flags_str",flagtext},
    {"check_outcomes","1=pass,0=fail,-1=not_evaluated"},
    {"se_note",mean?"Model-based standard error of the cosine ratio fit.":"Model-based standard errors on the log parameter scale."}};
  r.vectors.clear();for(auto &kv:vec)if(kv.first.rfind("windows.",0)==0)r.vectors[kv.first]=kv.second;
  if(!mean) {
    out["n_points"]=nan;
    for(auto k:{"t","t_scaled","cosine","cosine_se","cosine_fit"})r.setv(std::string("transform.")+k,{});
    if((method=="gaussian_jump" || method=="joint_voigt"))for(auto k:{"sine_scaled","sine_scaled_se","sine_fit"})r.setv(std::string("transform.")+k,{});
    double pilot=get("sigma_pilot");
    for(auto end:{"lo","hi"}) {auto v=r.vectors[std::string("windows.window_")+end+"_scaled"];for(double &x:v)x/=pilot;r.vectors[std::string("windows.window_")+end]=v;out[std::string("window_")+end]=get(std::string("u_")+end)/pilot;out[std::string("window_")+end+"_scaled"]=get(std::string("u_")+end);}
    if((method=="voigt" || method=="joint_voigt")) {
      r.vectors["windows.a"]=vec["cand_a_voigt"];r.vectors["windows.s2"]=vec["cand_s2_voigt"];
      r.vectors["windows.sigma_eq"]=r.vectors["windows.sigma"];
      auto v=vec["cand_s2_voigt"];for(auto &x:v)x=std::sqrt(x);r.vectors["windows.sigma"]=v;
      out["a"]=get("a_voigt");out["s2"]=get("s2_voigt");out.erase("a_voigt");out.erase("s2_voigt");
    }
    if((method=="gaussian_jump" || method=="joint_voigt"))r.vectors["windows.jump"]=vec["cand_aJ"];
  }
  r.scalars=std::move(out);
}
inline void add_simple(const sb_data* data,const sb_opts_mass &opts,sbc::Result &r) {
  sb_opts_mass o=opts;o.distribution_checks=0;o.check_all_windows=0;sb_result* h=nullptr;
  int rc=sb_mass_gaussian(data,&o,&h);
  if(rc)throw std::runtime_error("simple estimator calculation failed");
  auto *p=reinterpret_cast<sbc::Result*>(h);
  for(auto &kv:p->scalars)if(kv.first.rfind("loglinear_",0)==0)r.scalars[kv.first]=kv.second;
  sb_result_free(h);
}
inline int insufficient(sb_result** out,size_t n,const std::string &method,int threads=0) {
  auto *r=new sbc::Result;r->status=1;sbc::ThreadGuard tg(threads);r->set("threads_used",tg.count());r->set("n",n);r->set("n_input",n);r->set("flags",0);r->sets("reason","insufficient_data");
  for(auto k:{"B","sigma","mu","se","sigma_pilot","rank","tier","p","se_logB","se_logS"})r->set(k,NAN);
  finish(*r,method);r->finish_names();*out=reinterpret_cast<sb_result*>(r);return 0;
}
}
