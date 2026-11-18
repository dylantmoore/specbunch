#include "stplugin.h"
#include "specbunch.h"
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <memory>
#include <cstring>
#include <stdexcept>
#include <algorithm>
using Result = sb_result;
static std::unique_ptr<Result, decltype(&sb_result_free)> pending(nullptr, sb_result_free);
static char* mut(const std::string& x) { return const_cast<char*>(x.c_str()); }
static void check(int rc) { if(rc) throw std::runtime_error("Stata interface error " + std::to_string(rc)); }
static std::vector<std::string> names(const char* s) { std::vector<std::string> v; std::istringstream in(s?s:""); std::string k; while(in>>k) v.push_back(k); return v; }
static std::string alias(std::string s) { std::replace(s.begin(),s.end(),'.','_');  return s; }
static void local(std::string k,std::string v) { check(SF_macro_save(mut("_"+k),mut(v))); }
static std::string macro(std::string k) { std::vector<char> b(1000000); if(int rc=SF_macro_use(mut(k),b.data(),b.size()-1)) throw std::runtime_error("macro "+k+": "+std::to_string(rc)); return b.data(); }
static double real_argument(const char* text) {
 size_t used=0;double value=std::stod(text,&used);
 if(used!=std::strlen(text))throw std::runtime_error("invalid numeric plugin argument");
 return value;
}
static double stata(double d) { return std::isfinite(d) ? d : SV_missval; }
static std::string bits(double d) {
 uint64_t u; std::memcpy(&u,&d,8);
 std::string out(16, '0'); const char* hex="0123456789abcdef";
 for(int i=15;i>=0;--i) { out[i]=hex[u&15]; u >>= 4; }
 return out;
}
static void fail(std::string s) { s += "\n"; SF_error(mut(s)); }
static void publish(Result* r) {
 pending.reset(r);
  local("core_status",std::to_string(sb_result_status(r)));
  for(int kind=0;kind<3;++kind) { std::string list; int index=0; for(auto& k:names(sb_result_names(r,kind))) { alias(k); list+=k+" "; if(kind==1) { size_t n=0; sb_result_vector(r,k.c_str(),&n); local("len"+std::to_string(++index),std::to_string(n)); } if(kind==2) local("str"+std::to_string(++index),sb_result_string(r,k.c_str())); } local(kind==0?"scalar_names":kind==1?"vector_names":"string_names",list); }
}
STDLL
#if !defined(_WIN32)
__attribute__((visibility("default")))
#endif
stata_call(int argc, char* argv[]) {
 try {
  if(argc<1) throw std::runtime_error("Missing plugin operation");
  std::string op=argv[0];
  if(op=="load") { // Exact little-endian fixture transport, no estimation.
   if(argc!=SF_nvars()+1) throw std::runtime_error("load: column count mismatch");
   for(int j=1;j<argc;++j) { std::ifstream f(argv[j],std::ios::binary); if(!f) throw std::runtime_error("Cannot read fixture");
    for(int i=1;i<=SF_nobs();++i) { unsigned char b[8]; if(!f.read((char*)b,8)) throw std::runtime_error("Short fixture"); uint64_t u=0; for(int q=7;q>=0;--q) u=(u<<8)|b[q]; double d; std::memcpy(&d,&u,8); check(SF_vstore(j,i,stata(d))); }
    if(f.peek()!=EOF) throw std::runtime_error("Long fixture");
   } return 0;
  }
  if(op=="release") { pending.reset(); return 0; }
  if(op=="fetch_vectors") {
   if(!pending || argc!=2) throw std::runtime_error("fetch_vectors: missing result or matrix");
   auto keys=names(sb_result_names(pending.get(),1));
   if(SF_col(argv[1])!=(int)keys.size()) throw std::runtime_error("fetch_vectors: matrix column mismatch");
   int col=0;
   for(const auto& key:keys) {
    size_t n=0; const double* v=sb_result_vector(pending.get(),key.c_str(),&n);
    if(n>(size_t)SF_row(argv[1])) throw std::runtime_error("fetch_vectors: matrix row mismatch");
    ++col; std::string outbits; outbits.reserve(n*17);
    for(size_t row=0;row<n;++row) {
     check(SF_mat_store(argv[1],row+1,col,stata(v[row])));
     outbits+=bits(v[row]); outbits+=' ';
    }
    local("vbits"+std::to_string(col),outbits);
   }
   return 0;
  }
  if(op=="fetch") {
   if(!pending || argc!=3) throw std::runtime_error("fetch: missing result or arguments");
   std::string key=argv[1], outbits;
   if(key=="scalars") { int j=0; for(auto& k:names(sb_result_names(pending.get(),0))) { double d; check(sb_result_scalar(pending.get(),k.c_str(),&d)); check(SF_mat_store(argv[2],1,++j,stata(d))); outbits+=bits(d)+" "; } }
   else { size_t n=0; auto v=sb_result_vector(pending.get(),key.c_str(),&n); for(size_t j=0;j<n;++j) { check(SF_mat_store(argv[2],j+1,1,stata(v[j]))); outbits+=bits(v[j])+" "; } }
   local("bits",outbits); return 0;
  }
  if(op=="npmle") {
   if(argc!=18)throw std::runtime_error("invalid NPMLE plugin arguments");
   pending.reset();int nv=SF_nvars();if(nv<1||nv>2)throw std::runtime_error("NPMLE expects the running variable and optional weights");
   std::vector<double> z,w;for(int i=SF_in1();i<=SF_in2();++i)if(SF_ifobs(i)){double v;check(SF_vdata(1,i,&v));z.push_back(SF_is_missing(v)?NAN:v);if(nv==2){check(SF_vdata(2,i,&v));w.push_back(SF_is_missing(v)?NAN:v);}}
   sb_opts_npmle o;sb_opts_npmle_init(&o);o.B=real_argument(argv[1]);o.sigma=real_argument(argv[2]);o.a=real_argument(argv[3]);o.kink=real_argument(argv[4]);o.threads=std::stoi(argv[5]);
   o.core_intervals=std::stoi(argv[6]);o.grid_points=std::stoi(argv[7]);o.exact_threshold=std::stoi(argv[8]);o.max_iterations=std::stoi(argv[9]);
   o.lambda_min=real_argument(argv[10]);o.lambda_max=real_argument(argv[11]);o.lambda_points=std::stoi(argv[12]);o.refinement_points=std::stoi(argv[13]);
   o.histogram_width=real_argument(argv[14]);o.quadrature_tolerance=real_argument(argv[15]);o.optimality_tolerance=real_argument(argv[16]);
   std::vector<double> edges;if(std::string(argv[17])!="_none"){int nr=SF_row(argv[17]),nc=SF_col(argv[17]);if(nr!=1&&nc!=1)throw std::runtime_error("bins matrix must be a vector");for(int i=1;i<=nr;++i)for(int j=1;j<=nc;++j){double v;check(SF_mat_el(argv[17],i,j,&v));edges.push_back(v);}}
   sb_data d{z.size(),z.data(),nv==2?w.data():nullptr};Result* out=nullptr;int rc=sb_npmle(&d,&o,edges.empty()?nullptr:edges.data(),edges.size(),&out);
   if(rc){if(out)sb_result_free(out);throw std::runtime_error(sb_last_error());}publish(out);return 0;
  }
  if(op!="run" || argc!=13) throw std::runtime_error("run expects entry, center, threads, seed, policy, distribution_checks, candidates, budgets, band, p0, R");
  pending.reset();
  std::string entry=argv[1]; double center=real_argument(argv[2]); int threads=std::stoi(argv[3]);
  size_t used=0; uint64_t seed=std::stoull(argv[4],&used); if(used!=strlen(argv[4]) || argv[4][0]=='-') throw std::runtime_error("seed must be uint64 decimal");
  int base=false?3:(entry=="bmean"?2:1), nv=SF_nvars(); if(nv!=base && nv!=base+1) throw std::runtime_error("Wrong number of variables");
  std::vector<std::vector<double>> a(nv); for(auto& v:a) v.reserve(std::max(0,SF_in2()-SF_in1()+1));
  for(int i=SF_in1();i<=SF_in2();++i) if(SF_ifobs(i)) { for(int j=0;j<nv;++j) { double d; check(SF_vdata(j+1,i,&d)); if(SF_is_missing(d)) d=std::numeric_limits<double>::quiet_NaN(); if(j==0) d-=center; a[j].push_back(d); } }
  const double* wt=nv>base?a.back().data():nullptr; Result* r=nullptr; int rc;
if(entry=="bmean") { sb_opts_bmean o; sb_opts_bmean_init(&o); o.threads=threads; o.seed=seed; o.reserved_band=std::stoi(argv[10]); o.reserved_visibility=std::stoi(argv[11]); o.R=std::stoi(argv[12]); sb_data_xy d{a[0].size(),a[0].data(),a[1].data(),wt}; rc=sb_bmean(&d,&o,&r); }
  else { sb_opts_mass o; sb_opts_mass_init(&o); o.threads=threads; o.seed=seed; o.policy=std::stoi(argv[5]); o.distribution_checks=std::stoi(argv[6]); o.check_all_windows=std::stoi(argv[7]); o.lp_iterations=std::stoll(argv[8]); o.lp_retry_iterations=std::stoll(argv[9]); sb_data d{a[0].size(),a[0].data(),wt}; if(entry=="gaussian") rc=sb_mass_gaussian(&d,&o,&r); else if(entry=="voigt") rc=sb_mass_voigt(&d,&o,&r); else if(entry=="gaussian_jump") rc=sb_mass_gaussian_jump(&d,&o,&r); else if(entry=="joint_voigt") rc=sb_mass_joint_voigt(&d,&o,&r); else throw std::runtime_error("Unknown method"); }
  if(rc) { if(r) sb_result_free(r); fail("specbunch C error "+std::to_string(rc)+": "+sb_last_error()); return 498; }
  publish(r);
  return 0;
 } catch(const std::exception& e) { pending.reset(); fail(e.what()); return 498; } catch(...) { pending.reset(); fail("Unknown plugin failure"); return 498; }
}
