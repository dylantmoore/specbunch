#include "specbunch.h"
#include "sbc_internal.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

namespace sbpost {
constexpr double pi=3.14159265358979323846;
using Matrix=std::array<double,25>;
using Row=std::array<double,5>;
// Fixed-order Cholesky inverse for the small local regressions.
bool inverse(int k,Matrix& A){
 Matrix L{},Li{};
 for(int i=0;i<k;++i)for(int j=0;j<=i;++j){double s=A[i*5+j];for(int p=0;p<j;++p)s-=L[i*5+p]*L[j*5+p];if(i==j){if(!(s>0))return false;L[i*5+i]=std::sqrt(s);}else L[i*5+j]=s/L[j*5+j];}
 for(int i=0;i<k;++i){Li[i*5+i]=1/L[i*5+i];for(int j=0;j<i;++j){double s=0;for(int p=j;p<i;++p)s+=L[i*5+p]*Li[p*5+j];Li[i*5+j]=-s/L[i*5+i];}}
 for(int i=0;i<k;++i)for(int j=0;j<k;++j){double s=0;for(int p=std::max(i,j);p<k;++p)s+=Li[p*5+i]*Li[p*5+j];A[i*5+j]=s;}return true;
}
struct Law {
 double sigma,a;mutable bool uncertain=false;mutable double error=0;
 double cdf(double x)const{
  if(a==0)return .5*std::erfc(-x/(sigma*std::sqrt(2.)));
  if(sigma==0)return .5+std::atan(x/a)/pi;
  auto integrate=[&](int n){double L=10*sigma,h=2*L/n,sum=0;for(int i=0;i<=n;++i){double g=-L+i*h,w=(i==0||i==n)?1:(i%2?4:2);sum+=w*std::exp(-.5*g*g/(sigma*sigma))*(.5+std::atan((x-g)/a)/pi);}return sum*h/3/(sigma*std::sqrt(2*pi));};
  double old=integrate(2000),value=old;
  for(int n=4000;n<=64000;n*=2){value=integrate(n);double e=std::fabs(value-old)/15;error=std::max(error,e);if(e<=1e-11)return value;old=value;}
  uncertain=true;return value;
 }
 double cf(double t)const{return std::exp(-.5*sigma*sigma*t*t-a*std::fabs(t));}
};
bool coherent(double B,double s,double a,double kink){return B>0&&B<1&&s>=0&&a>=0&&s+a>0&&std::isfinite(B+s+a+kink);}
void flag(sbc::Result& r,const char* name,bool value,std::string& flags){r.set(name,value);if(value){if(!flags.empty())flags+=' ';flags+=name;}}
int finish(std::unique_ptr<sbc::Result>& r,sb_result** out,const char* reason=""){
 if(*reason && r->scalars.count("saez.e")){r->sets("saez.status",reason);r->sets("continuous.status",reason);}
 r->sets("refusal_reason",reason);r->sets("status",*reason?"refused":"ok");r->status=*reason?1:0;r->set("status_code",r->status);r->finish_names();*out=reinterpret_cast<sb_result*>(r.release());return SBC_OK;
}
struct Bin {double w=0,w2=0,wx=0,mean=NAN,variance=NAN,center=0,atom=0,share=NAN;};
int nonbuncher(const sb_data_xy* d,const sb_opts_nonbuncher_mean* supplied,sb_result** out){
 if(!out)return sbc::fail(SBC_E_NULL,"result pointer is NULL");*out=nullptr;
 if(!d||!d->z||!d->y||!supplied)return sbc::fail(SBC_E_NULL,"sample, characteristic and options are required");
 if(int rc=sbc::check_opts_header(supplied,8,SB_ABI_VERSION))return rc;
 sb_opts_nonbuncher_mean o;sb_opts_nonbuncher_mean_init(&o);std::memcpy(&o,supplied,std::min<size_t>(sizeof o,supplied->struct_size));
 if(!coherent(o.B,o.sigma,o.a,o.kink)||!std::isfinite(o.mu_b)||!(o.se_mu_b>=0)||!std::isfinite(o.se_mu_b)||!(o.bandwidth>0)||!(o.bin_width>0)||!std::isfinite(o.bandwidth+o.bin_width)||(o.degree!=1&&o.degree!=2))return sbc::fail(SBC_E_INPUT,"invalid mass, friction, buncher mean or local polynomial options");
 double H=std::max(o.bandwidth,.10),extent=std::ceil(H/o.bin_width);
 if(!(extent>=2&&extent<=10000))return sbc::fail(SBC_E_INPUT,"bin_width must yield between 4 and 20000 smoothing bins");
 int half=(int)extent;std::vector<Bin> bins(2*half);double W=0,W2=0,xmin=INFINITY,xmax=-INFINITY;size_t count=0;
 sbc::ThreadGuard threads(o.threads);
 auto usable=[&](size_t i){return std::isfinite(d->z[i]-o.kink)&&std::isfinite(d->y[i])&&(!d->w||(std::isfinite(d->w[i])&&d->w[i]>0));};
 for(size_t i=0;i<d->n;++i)if(usable(i)){double z=d->z[i]-o.kink,w=d->w?d->w[i]:1,x=d->y[i];W+=w;W2+=w*w;++count;xmin=std::min(xmin,x);xmax=std::max(xmax,x);double j=std::floor(z/o.bin_width);if(j>=-half&&j<half){auto& b=bins[(int)j+half];b.w+=w;b.w2+=w*w;b.wx+=w*x;}}
 auto result=std::make_unique<sbc::Result>();auto& r=*result;
 r.set("n",count);r.set("n_input",d->n);r.set("n_eff",W2>0?W*W/W2:0);r.set("threads_used",threads.count());r.set("bandwidth",o.bandwidth);r.set("bin_width",o.bin_width);r.set("degree",o.degree);r.set("kink",o.kink);r.set("B",o.B);r.set("sigma",o.sigma);r.set("a",o.a);r.set("mu_b",o.mu_b);r.set("se_mu_b",o.se_mu_b);r.set("m_NB",NAN);r.set("se",NAN);r.set("pi_kink",NAN);r.sets("flags","");
 r.sets("estimator","triangular-kernel local polynomial with fixed buncher mean");
 r.sets("se_note","Standard error, not a confidence interval; conditional on mass, friction and smoothed density; buncher-mean variance propagated assuming independence.");
 if(count<10||!(W>0&&std::isfinite(W+W2)))return finish(result,out,"insufficient_data");
 for(auto& b:bins)if(b.w>0){b.mean=b.wx/b.w;b.variance=0;}
 for(size_t i=0;i<d->n;++i)if(usable(i)){double j=std::floor((d->z[i]-o.kink)/o.bin_width);if(j>=-half&&j<half){auto& b=bins[(int)j+half];double e=d->y[i]-b.mean;b.variance+=(d->w?d->w[i]:1)*e*e;}}
 // Kish-corrected weighted within-bin variance of the mean. The density
 // denominator uses every usable observation, including those outside H.
 Law law{o.sigma,o.a};Matrix A{};Row rhs{};std::vector<Row> design(bins.size());int density_bins=0;
 for(size_t j=0;j<bins.size();++j){auto& b=bins[j];double lo=((int)j-half)*o.bin_width,hi=lo+o.bin_width;b.center=(lo+hi)/2;
  if(!(b.w>0)||std::fabs(b.center)>H)continue;
  double den=b.w-b.w2/b.w;b.variance=den>0?b.variance/den*b.w2/(b.w*b.w):0;
  b.atom=o.B*(law.cdf(hi)-law.cdf(lo))/o.bin_width;
  design[j]={1,b.center,b.center*b.center,law.cdf(b.center)-.5,0};double y=b.w/(W*o.bin_width)-b.atom;++density_bins;
  for(int p=0;p<4;++p){rhs[p]+=design[j][p]*y;for(int q=0;q<4;++q)A[p*5+q]+=design[j][p]*design[j][q];}
 }
 if(density_bins<6||!inverse(4,A))return finish(result,out,"insufficient_local_variation");
 Row density{};for(int p=0;p<4;++p)for(int q=0;q<4;++q)density[p]+=A[p*5+q]*rhs[q];
 bool clipped=false,singleton=false;std::vector<Row> X;std::vector<double> weight,y,var,share,centers,means,binweights;
 for(size_t j=0;j<bins.size();++j){auto& b=bins[j];if(!(b.w>0)||std::fabs(b.center)>=o.bandwidth)continue;double f=0;for(int p=0;p<4;++p)f+=density[p]*design[j][p];clipped|=f<=1e-12;b.share=b.atom/(b.atom+std::max(f,1e-12));double q=1-b.share;
  X.push_back({b.share,q,q*b.center,o.degree==2?q*b.center*b.center:0,0});weight.push_back(b.w*(1-std::fabs(b.center)/o.bandwidth));y.push_back(b.mean-o.mu_b*b.share);var.push_back(b.variance);share.push_back(b.share);centers.push_back(b.center+o.kink);means.push_back(b.mean);binweights.push_back(b.w);singleton|=b.w*b.w<=b.w2;
 }
 int k=o.degree+2;r.set("n_bins",X.size());
 if(X.size()<(size_t)(k+2))return finish(result,out,"insufficient_data");
 A={};rhs={};for(size_t j=0;j<X.size();++j)for(int p=0;p<k;++p){rhs[p]+=weight[j]*X[j][p]*y[j];for(int q=0;q<k;++q)A[p*5+q]+=weight[j]*X[j][p]*X[j][q];}
 for(int p=0;p<k;++p){A[p*5]=0;A[p]=0;}A[0]=1;rhs[0]=0;
 if(!inverse(k,A))return finish(result,out,"insufficient_local_variation");
 Row beta{},derivative{},t{};for(int p=0;p<k;++p)for(int q=0;q<k;++q)beta[p]+=A[p*5+q]*rhs[q];
 Matrix meat{},V{};for(size_t j=0;j<X.size();++j){double c=weight[j]*weight[j]*var[j];for(int p=1;p<k;++p){t[p]+=weight[j]*X[j][p]*share[j];for(int q=1;q<k;++q)meat[p*5+q]+=c*X[j][p]*X[j][q];}}
 for(int p=0;p<k;++p)for(int q=0;q<k;++q)for(int u=0;u<k;++u)for(int v=0;v<k;++v)V[p*5+q]+=A[p*5+u]*meat[u*5+v]*A[v*5+q];
 for(int p=1;p<k;++p)for(int q=1;q<k;++q)derivative[p]-=A[p*5+q]*t[q];
 double se=std::sqrt(std::max(0.,V[6]+derivative[1]*derivative[1]*o.se_mu_b*o.se_mu_b));
 // Evaluate the same bin-averaged density model in a bin centred at the kink.
 double atom0=o.B*(law.cdf(o.bin_width/2)-law.cdf(-o.bin_width/2))/o.bin_width;
 r.set("pi_kink",atom0/(atom0+std::max(density[0],1e-12)));r.set("m_NB",beta[1]);r.set("se",se);r.set("density_bandwidth",H);r.set("quadrature_error",law.error);
 std::vector<double> fitted,nb;for(size_t j=0;j<X.size();++j){double v=0;for(int p=1;p<k;++p)v+=X[j][p]*beta[p];fitted.push_back(o.mu_b*share[j]+v);double z=centers[j]-o.kink;nb.push_back(beta[1]+beta[2]*z+(o.degree==2?beta[3]*z*z:0));}
 r.setv("local_fit.center",centers);r.setv("local_fit.observed_mean",means);r.setv("local_fit.fitted_mixture",fitted);r.setv("local_fit.nonbuncher_mean",nb);r.setv("local_fit.buncher_share",share);r.setv("local_fit.weight",binweights);
 std::string flags;flag(r,"outside_outcome_range",beta[1]<xmin||beta[1]>xmax,flags);flag(r,"density_floor_used",clipped||density[0]<=1e-12,flags);flag(r,"sparse_bins",singleton,flags);flag(r,"quadrature_inconclusive",law.uncertain,flags);r.sets("flags",flags);
 if(!std::isfinite(beta[1]+se))return finish(result,out,"numerical_failure");return finish(result,out);
}
// A fixed grid and fixed chunk reductions make the weighted empirical CF
// independent of scheduling and thread count. No platform-dependent long doubles.
void ecf(const std::vector<double>& z,const std::vector<double>& w,double delta,int M,std::vector<double>& C,std::vector<double>& S){
 constexpr int block=512,chunks=64;size_t n=z.size(),nb=(n+block-1)/block,per=(nb+chunks-1)/chunks;
 std::vector<double> partial((size_t)chunks*2*M,0),totals(chunks,0);
 #pragma omp parallel for schedule(static)
 for(int ch=0;ch<chunks;++ch){std::array<double,block> c{},s{},cd{},sd{};double* p=&partial[(size_t)ch*2*M];
  for(size_t b=ch*per;b<std::min(nb,(ch+1)*per);++b){size_t start=b*block,len=std::min((size_t)block,n-start);for(size_t i=0;i<len;++i){c[i]=1;s[i]=0;cd[i]=std::cos(delta*z[start+i]);sd[i]=std::sin(delta*z[start+i]);totals[ch]+=w[start+i];}
   for(int j=0;j<M;++j){if(j)for(size_t i=0;i<len;++i){double old=c[i];c[i]=old*cd[i]-s[i]*sd[i];s[i]=s[i]*cd[i]+old*sd[i];}double qc=0,qs=0;for(size_t i=0;i<len;++i){qc+=w[start+i]*c[i];qs+=w[start+i]*s[i];}p[j]+=qc;p[M+j]+=qs;}
  }
 }
 C.assign(M,0);S.assign(M,0);double W=0;for(int ch=0;ch<chunks;++ch){W+=totals[ch];for(int j=0;j<M;++j){C[j]+=partial[(size_t)ch*2*M+j];S[j]+=partial[(size_t)ch*2*M+M+j];}}for(int j=0;j<M;++j){C[j]/=W;S[j]/=W;}
}
int elasticity(const sb_data* d,const sb_opts_eti* supplied,sb_result** out){
 if(!out)return sbc::fail(SBC_E_NULL,"result pointer is NULL");*out=nullptr;if(!d||!d->z||!supplied)return sbc::fail(SBC_E_NULL,"sample and options are required");
 if(int rc=sbc::check_opts_header(supplied,8,SB_ABI_VERSION))return rc;
 sb_opts_eti o;sb_opts_eti_init(&o);std::memcpy(&o,supplied,std::min<size_t>(sizeof o,supplied->struct_size));
 if(!coherent(o.B,o.sigma,o.a,o.kink)||!(o.window_lo>0&&o.window_hi>o.window_lo)||!std::isfinite(o.window_hi)||!(o.net_below>o.net_above&&o.net_above>0)||!std::isfinite(o.net_below)||o.integration_intervals<4||o.integration_intervals>4000)return sbc::fail(SBC_E_INPUT,"ETI requires coherent mass/friction, a positive frequency band and net_below > net_above > 0");
 int intervals=o.integration_intervals+(o.integration_intervals%2);double delta=o.window_lo/intervals,points=std::ceil(o.window_hi/delta)+1;
 if(!(points<=20001))return sbc::fail(SBC_E_INPUT,"ETI frequency grid exceeds 20001 points");int M=(int)points;
 sbc::ThreadGuard threads(o.threads);std::vector<double> z,w;double W=0,W2=0;
 for(size_t i=0;i<d->n;++i){double zi=d->z[i]-o.kink,wi=d->w?d->w[i]:1;if(!std::isfinite(zi)||!std::isfinite(wi)||!(wi>0))continue;z.push_back(zi);w.push_back(wi);W+=wi;W2+=wi*wi;}
 auto result=std::make_unique<sbc::Result>();auto& r=*result;double q=std::log(o.net_below/o.net_above);
 r.set("n",z.size());r.set("n_input",d->n);r.set("n_eff",W2>0?W*W/W2:0);r.set("threads_used",threads.count());r.set("B",o.B);r.set("sigma",o.sigma);r.set("a",o.a);r.set("kink",o.kink);r.set("window_lo",o.window_lo);r.set("window_hi",o.window_hi);r.set("net_below",o.net_below);r.set("net_above",o.net_above);r.set("log_net_ratio",q);r.set("integration_intervals",intervals);r.sets("flags","");r.sets("estimator","shift-and-complete with empirical characteristic-function shoulders");
 for(auto name:{"saez.e","saez.d","continuous.e","continuous.d","shoulder_left","shoulder_right","shoulder_midpoint","shoulder_jump"})r.set(name,NAN);
 if(z.size()<2||!std::isfinite(W+W2)||!(W>0))return finish(result,out,"insufficient_data");
 Law law{o.sigma,o.a};if(!(law.cf(o.window_hi)>1e-150))return finish(result,out,"friction_underflow");
 std::vector<double> C,S;ecf(z,w,delta,M,C,S);double integral=0;
 for(int j=0;j<=intervals;++j){double t=j*delta,weight=(j==0||j==intervals)?1:(j%2?4:2);integral+=weight*C[j]/law.cf(t);}
 double midpoint=(integral*delta/3-o.B*o.window_lo)/pi,sw=0,swy=0;
 for(int j=intervals;j<M;++j){double t=j*delta;if(t>o.window_hi*(1+1e-12))break;double ph=law.cf(t),y=t*S[j]/ph,weight=ph*ph/(t*t);sw+=weight;swy+=weight*y;}
 double jump=swy/sw,left=midpoint-jump/2,right=midpoint+jump/2;
 r.set("shoulder_midpoint",midpoint);r.set("shoulder_jump",jump);r.set("shoulder_left",left);r.set("shoulder_right",right);
 std::string flags;flag(r,"negative_shoulder",left<0||right<0,flags);flag(r,"zero_shoulder",left==0||right==0,flags);r.sets("flags",flags);
 if(!(midpoint>0)||!std::isfinite(midpoint+jump))return finish(result,out,"nonpositive_density");
 double ds=o.B/midpoint;r.set("saez.d",ds);r.set("saez.e",ds/q);r.sets("saez.status","ok");
 if(left>=0&&right>=0){double dc=3*o.B/(left+right+std::sqrt(left*right));r.set("continuous.d",dc);r.set("continuous.e",dc/q);r.sets("continuous.status","ok");}else r.sets("continuous.status","negative_shoulder");
 return finish(result,out);
}
}
extern "C" {
SB_API void sb_opts_nonbuncher_mean_init(sb_opts_nonbuncher_mean* o){if(!o)return;std::memset(o,0,sizeof *o);o->struct_size=sizeof *o;o->abi_version=SB_ABI_VERSION;o->degree=1;o->B=NAN;o->sigma=NAN;o->mu_b=NAN;o->bandwidth=.05;o->bin_width=.0025;}
SB_API void sb_opts_eti_init(sb_opts_eti* o){if(!o)return;std::memset(o,0,sizeof *o);o->struct_size=sizeof *o;o->abi_version=SB_ABI_VERSION;o->integration_intervals=200;o->B=NAN;o->sigma=NAN;o->window_lo=NAN;o->window_hi=NAN;o->net_below=NAN;o->net_above=NAN;}
SB_API int sb_nonbuncher_mean(const sb_data_xy* d,const sb_opts_nonbuncher_mean* o,sb_result** out){return sbc::guarded([&]{return sbpost::nonbuncher(d,o,out);});}
SB_API int sb_eti(const sb_data* d,const sb_opts_eti* o,sb_result** out){return sbc::guarded([&]{return sbpost::elasticity(d,o,out);});}
}
