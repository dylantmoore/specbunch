#include "transform_output.hpp"
#include "public_result.hpp"








#include "llcore_prelude.hpp"  
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

namespace sbjoint_reference {
#define main sbjoint_reference_main_unused
#include "cfj/cfj.cpp"
#undef main
}  

#include "policy_run.hpp"
#include "joint_voigt_engine.hpp"

namespace {
int mass_joint_voigt_impl(const sb_data* data,const sb_opts_mass* supplied,sb_result** out) {
  if(!data || !out)return sbc::fail(SBC_E_NULL,"data or out is NULL");
  *out=nullptr;sb_opts_mass o;if(int rc=sbrun::read_opts(supplied,o))return rc;
  if(data->w)return sbc::fail(SBC_E_INPUT,"joint does not support weights");
  if(data->n==0)return sbpublic::insufficient(out,0,"joint_voigt",o.threads);
  if(!data->z)return sbc::fail(SBC_E_INPUT,"sample pointer is NULL");
  std::vector<double> z;z.reserve(data->n);
  for(size_t i=0;i<data->n;++i)if(std::isfinite(data->z[i]))z.push_back(data->z[i]);
  if(z.size()<100)return sbpublic::insufficient(out,z.size(),"joint_voigt",o.threads);
  auto result=std::make_unique<sbc::Result>();auto* r=result.get();sbjointvoigt::Output v;int selected=-1;
  {
    sbc::ThreadGuard tg(o.threads);r->set("threads_used",tg.count());v=sbjointvoigt::run(z,o.seed);
    auto frequency=[&](int j){return j>=0 && j<(int)v.u.size()?v.u[j]:NAN;};
    sb_opts_mass checks=o;checks.distribution_checks=0;
    selected=sbrun::run_policy(r,v.table,10,v.reason,nullptr,checks,frequency,"joint_voigt",false);
  }
  r->set("sigma_eq",r->scalars["sigma"]);r->set("se_log_sigma_eq",r->scalars["se_logS"]);
  r->set("a_voigt",selected>=0?v.parameters[(size_t)selected*2]:NAN);
  r->set("s2_voigt",selected>=0?v.parameters[(size_t)selected*2+1]:NAN);
  r->set("sigma",std::sqrt(r->scalars["s2_voigt"]));
  r->set("aJ",selected>=0?v.table[(size_t)selected*10+9]:NAN);
  r->set("n",z.size());r->set("n_input",data->n);r->set("weighted",0);r->set("sigma_pilot",v.pilot);r->set("seed",(double)o.seed);
  size_t count=v.table.size()/10;std::vector<double> a(count),s2(count),jump(count);
  for(size_t i=0;i<count;++i){a[i]=v.parameters[2*i];s2[i]=v.parameters[2*i+1];jump[i]=v.table[10*i+9];}
  r->setv("cand_a_voigt",std::move(a));r->setv("cand_s2_voigt",std::move(s2));r->setv("cand_aJ",std::move(jump));
  sbpublic::add_simple(data,o,*r);sbpublic::finish(*r,"joint_voigt");sbtransform::add(*data,o.threads,*r);r->finish_names();
  *out=reinterpret_cast<sb_result*>(result.release());return SBC_OK;
}
}
extern "C" SB_API int sb_mass_joint_voigt(const sb_data* data,const sb_opts_mass* opts,sb_result** out) {
  return sbc::guarded([&]{return mass_joint_voigt_impl(data,opts,out);});
}
