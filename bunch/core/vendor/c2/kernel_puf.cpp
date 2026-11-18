// kernel_puf.cpp: the frozen C2 screen kernel (cjms-estimators/cppest/refuse/src_frozen/full/kernel.cpp, sha256 d20767bd1c23f786...) for the PUF run, 2026-11-06.
// One change: moment_cuts (LAPACK zheevd; used only by solver.Exchange, which the C2 screen never calls) is compiled only
// when P1_WITH_LAPACK is defined, else it returns -999, so the library builds on Linux without LAPACK. cdf_screen,
// moment_box and box_support (the functions C2 uses) are unchanged.
#ifdef P1_WITH_LAPACK
#include <Accelerate/Accelerate.h>
#endif
#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>
// Column-major Hermitian Toeplitz, upper diagonal q_h=c_h+i*s_h.
#ifdef P1_WITH_LAPACK
extern "C" int moment_cuts(int K,const double* q,double mass,int want,double* eig,double* cuts,double* vectors){
 int N=K+1,info=0,lda=N; char job='V',uplo='U';
 std::vector<__CLPK_doublecomplex>A(N*N);
 for(int j=0;j<N;++j)for(int i=0;i<N;++i){int h=std::abs(j-i);A[i+j*N].r=h?q[h-1]:mass;A[i+j*N].i=h?(j>i?q[K+h-1]:-q[K+h-1]):0;}
 int lw=-1,lr=-1,li=-1,iq;__CLPK_doublecomplex wq;double rq;
 zheevd_(&job,&uplo,&N,A.data(),&lda,eig,&wq,&lw,&rq,&lr,&iq,&li,&info);
 if(info)return info;lw=int(wq.r);lr=int(rq);li=iq;
 std::vector<__CLPK_doublecomplex>w(lw);std::vector<double>rw(lr);std::vector<int>iw(li);
 zheevd_(&job,&uplo,&N,A.data(),&lda,eig,w.data(),&lw,rw.data(),&lr,iw.data(),&li,&info);
 if(info)return info;
 for(int k=0;k<std::min(want,N);++k){
  for(int i=0;i<N;++i){vectors[k*2*N+i]=A[i+k*N].r;vectors[k*2*N+N+i]=A[i+k*N].i;}
  for(int h=1;h<=K;++h){double re=0,im=0;for(int i=0;i<N-h;++i){auto a=A[i+k*N],b=A[i+h+k*N];re+=a.r*b.r+a.i*b.i;im+=a.r*b.i-a.i*b.r;}
   cuts[k*2*K+h-1]=2*re;cuts[k*2*K+K+h-1]=-2*im;}
 }
 return 0;
}
#else
extern "C" int moment_cuts(int, const double*, double, int, double*, double*, double*) { return -999; }
#endif
// Box for unnormalized background moments, without amplifying small d.
extern "C" int moment_box(int K,const double*c,const double*s,int n,double B,double sigma,double*lo,double*hi,double*diskmargin){
 double r=std::sqrt(2*std::log(512/.01)/n),m=1-B,margin=-1;
 for(int h=0;h<K;++h){double t=(h+1)*.25,d=std::exp(-.5*(sigma*t)*(sigma*t));
 double cl=std::max(-1.,c[h]-r),cu=std::min(1.,c[h]+r),sl=std::max(-1.,s[h]-r),su=std::min(1.,s[h]+r);
 double dx=std::max({cl-B*d,B*d-cu,0.}),dy=std::max({sl,-su,0.});margin=std::max(margin,std::hypot(dx,dy)-m*d);
 // Clip before division; if d=0, coordinate exclusion detects nonzero bands.
 auto lower=[&](double x,double shift){return d==0?-m:(x<=d*(shift-m)?-m:(x>=d*(shift+m)?m:x/d-shift));};
 auto upper=[&](double x,double shift){return d==0?m:(x>=d*(shift+m)?m:(x<=d*(shift-m)?-m:x/d-shift));};
 lo[h]=lower(cl,B);hi[h]=upper(cu,B);lo[K+h]=lower(sl,0);hi[K+h]=upper(su,0);
 }
 *diskmargin=margin;return margin>1e-12?1:0;
}
extern "C" double box_support(int d,const double*a,const double*lo,const double*hi,double mass){
 double v=mass;for(int i=0;i<d;++i)v+=a[i]>=0?a[i]*hi[i]:a[i]*lo[i];return v;
}
extern "C" int cdf_screen(const double*x,const double*f,int n,double B,double sig,double*margin,int*which){
 constexpr int J=66;double eps=std::sqrt(std::log(2/.001)/(2*n)),best=-1;int hit=-1;double G[J];
 for(int j=0;j<J;++j){G[j]=sig>0?.5*std::erfc(-x[j]/(sig*std::sqrt(2.))):(x[j]>=0?1.:0.);
  double a=B*G[j]-(f[j]+eps),b=(f[j]-eps)-(B*G[j]+1-B);if(std::max(a,b)>best){best=std::max(a,b);hit=j;}}
 for(int j=0;j<J;++j)for(int i=0;i<j;++i){double atom=B*(G[j]-G[i]),width=x[j]-x[i];double maxp=sig>0?std::erf(width/(2*sig*std::sqrt(2.))):1.;double upper=atom+(1-B)*maxp;
  double v=std::max((f[j]-f[i]-2*eps)-upper,atom-(f[j]-f[i]+2*eps));if(v>best){best=v;hit=J+i*J+j;}}
 *margin=best;*which=hit;return best>1e-10;
}
