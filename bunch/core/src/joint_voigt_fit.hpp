#pragma once
#include <array>
#include <cmath>
#include <algorithm>
#include <limits>
namespace sbjointvoigt {
// Work in u=t*sigma_pilot. alpha=a/sigma_pilot, beta=s2/sigma_pilot^2,
// J=jump*sigma_pilot. Symmetric additive friction convolves the entire latent
// distribution: its real characteristic function multiplies BOTH cosine and
// sine. Thus the density-jump term is J*exp(-alpha*u-beta*u*u/2)/u.
struct Fit {
  bool valid=false;
  double B=NAN,J=NAN,alpha=NAN,beta=NAN,objective=INFINITY,D=NAN,seB=NAN,seEq=NAN;
};
struct Profile {
  bool valid=false; double B=0,J=0,f=0,g[4]={},H[4][4]={};
};
inline bool solve(double H[4][4],const double* g,int n,double* x) {
  double A[4][5]={};
  for(int i=0;i<n;++i){for(int j=0;j<n;++j)A[i][j]=H[i][j];A[i][n]=g[i];}
  for(int j=0;j<n;++j) {
    int pivot=j;for(int i=j+1;i<n;++i)if(std::abs(A[i][j])>std::abs(A[pivot][j]))pivot=i;
    if(!(std::abs(A[pivot][j])>1e-24))return false;
    for(int k=j;k<=n;++k)std::swap(A[j][k],A[pivot][k]);
    double v=A[j][j];for(int k=j;k<=n;++k)A[j][k]/=v;
    for(int i=0;i<n;++i)if(i!=j){v=A[i][j];for(int k=j;k<=n;++k)A[i][k]-=v*A[j][k];}
  }
  for(int i=0;i<n;++i)x[i]=A[i][n];return true;
}
template<class Context,class Moments>
inline Profile profile(const Context& c,const Moments& M,const double* C,const double* S,int first,int count,double alpha,double beta) {
  Profile p;double h00=0,h01=0,h11=0,b0=0,b1=0,E[160],u[160];
  const double invn=1./c.n;
  for(int k=0;k<count;++k) {
    int j=first+k;u[k]=c.t[j]*c.sigma_p;E[k]=std::exp(-alpha*u[k]-.5*beta*u[k]*u[k]);
    double e=E[k],v=e/u[k],kap=M.kap[j],w1=M.w1[j]*invn,w2=M.w2[j]*invn,y=C[j]-kap*S[j];
    h00+=w1*e*e;h01-=w1*e*kap*v;h11+=(w1*kap*kap+w2)*v*v;
    b0+=w1*e*y;b1+=-w1*kap*v*y+w2*v*S[j];
  }
  double det=h00*h11-h01*h01;
  if(!(det>1e-24*std::max(h00*h11,1e-300)))return p;
  p.B=(b0*h11-b1*h01)/det;p.J=(b1*h00-b0*h01)/det;
  if(!std::isfinite(p.B)||!std::isfinite(p.J))return p;
  for(int k=0;k<count;++k) {
    int j=first+k;double e=E[k],v=e/u[k],kap=M.kap[j],w1=M.w1[j]*invn,w2=M.w2[j]*invn;
    double m1=p.B*e-kap*p.J*v,m2=p.J*v,r1=C[j]-kap*S[j]-m1,r2=S[j]-m2;
    double d1[4]={e,-kap*v,-u[k]*m1,-.5*u[k]*u[k]*m1};
    double d2[4]={0,v,-u[k]*m2,-.5*u[k]*u[k]*m2};
    p.f+=w1*r1*r1+w2*r2*r2;
    for(int i=0;i<4;++i){p.g[i]+=w1*r1*d1[i]+w2*r2*d2[i];for(int h=0;h<4;++h)p.H[i][h]+=w1*d1[i]*d1[h]+w2*d2[i]*d2[h];}
  }
  p.valid=std::isfinite(p.f);return p;
}
template<class Context,class Moments>
inline Fit optimize(const Context& c,const Moments& M,const double* C,const double* S,int first,int count,double alpha,double beta,int boundary) {
  Fit out;if(boundary==1)alpha=0;if(boundary==2)beta=0;
  Profile p=profile(c,M,C,S,first,count,alpha,beta);
  if(!p.valid)return out;
  bool converged=false;
  for(int it=0;it<100;++it) {
    const double det=p.H[0][0]*p.H[1][1]-p.H[0][1]*p.H[0][1];
    if(!(det>0))break;
    double H[2][2];
    for(int i=0;i<2;++i)for(int j=0;j<2;++j) {
      double x=(p.H[1][1]*p.H[0][j+2]-p.H[0][1]*p.H[1][j+2])/det;
      double y=(p.H[0][0]*p.H[1][j+2]-p.H[0][1]*p.H[0][j+2])/det;
      H[i][j]=p.H[i+2][j+2]-p.H[i+2][0]*x-p.H[i+2][1]*y;
    }
    double g[2]={p.g[2],p.g[3]},x[2]={alpha,beta},scaled=0;
    bool active[2]={boundary==1 || (alpha==0 && g[0]<=0),boundary==2 || (beta==0 && g[1]<=0)};
    for(int j=0;j<2;++j)if(!active[j])scaled=std::max(scaled,std::abs(g[j])/std::sqrt(std::max(p.H[j+2][j+2],1e-24)));
    if(scaled<=1e-10*(1+std::sqrt(p.f))){converged=true;break;}
    double step[2]={};bool accepted=false;
    for(int ridge=0;ridge<10 && !accepted;++ridge) {
      double h00=std::max(H[0][0],0.)+std::pow(10.,ridge-10)*(1+p.H[2][2]);
      double h11=std::max(H[1][1],0.)+std::pow(10.,ridge-10)*(1+p.H[3][3]);
      double determinant=h00*h11-H[0][1]*H[1][0];
      if(active[0])step[1]=active[1]?0:g[1]/h11;
      else if(active[1])step[0]=g[0]/h00;
      else if(determinant>0){step[0]=(h11*g[0]-H[0][1]*g[1])/determinant;step[1]=(h00*g[1]-H[1][0]*g[0])/determinant;}
      else continue;
      double fraction=std::min(1.,2./std::max({std::abs(step[0]),std::abs(step[1]),1e-300}));
      for(int ls=0;ls<30;++ls,fraction*=.5) {
        double a=std::max(0.,x[0]+fraction*step[0]),b=std::max(0.,x[1]+fraction*step[1]);
        Profile q=profile(c,M,C,S,first,count,a,b);
        double next_scaled=0.;
        for(int j=0;j<2;++j)if(!((j==0?boundary==1:boundary==2) || ((j==0?a:b)==0 && q.g[j+2]<=0)))
          next_scaled=std::max(next_scaled,std::abs(q.g[j+2])/std::sqrt(std::max(q.H[j+2][j+2],1e-24)));
        const double roundoff=8*std::numeric_limits<double>::epsilon()*std::max(1.,p.f);
        if(q.valid && (q.f<p.f || (q.f<=p.f+roundoff && next_scaled<scaled))) {alpha=a;beta=b;p=q;accepted=true;break;}
      }
    }
    if(!accepted) {
      // Near the optimum an objective decrease can be below one rounding unit.
      converged=scaled<=1e-8*(1+std::sqrt(p.f));break;
    }
  }
  if(!converged)return out;
  for(int j=0;j<2;++j)if((j==0?alpha:beta)==0 && p.g[j+2]>1e-8*(1+std::sqrt(p.f))*std::sqrt(std::max(p.H[j+2][j+2],1e-24)))return out;
  out.valid=true;out.B=p.B;out.J=p.J;out.alpha=alpha;out.beta=beta;out.objective=p.f;return out;
}
template<class Context,class Moments>
inline Fit fit(const Context& c,const Moments& M,const double* C,const double* S,int first,int count,const Fit* initial,bool want_se) {
  Fit best;
  const double a=initial?initial->alpha:.2,b=initial?initial->beta:.8;
  for(int boundary=0;boundary<3;++boundary) {
    Fit q=optimize(c,M,C,S,first,count,a,b,boundary);
    if(q.valid && q.objective<best.objective)best=q;
    if(initial && q.valid)break;
  }
  if(!best.valid && !initial) {
    for(auto start:{std::array<double,2>{0.,1.},std::array<double,2>{1.,0.}}) {
      Fit q=optimize(c,M,C,S,first,count,start[0],start[1],0);
      if(q.valid && q.objective<best.objective)best=q;
    }
  }
  if(!best.valid)return best;
  best.D=0;
  for(int j=first;j<first+count;++j) {
    double u=c.t[j]*c.sigma_p,e=std::exp(-best.alpha*u-.5*best.beta*u*u);
    best.D=std::max(best.D,std::abs(C[j]-best.B*e)/c.s[j]);
    best.D=std::max(best.D,std::abs(S[j]-best.J*e/u)/M.sS[j]);
  }
  if(want_se && best.B>0) {
    Profile p=profile(c,M,C,S,first,count,best.alpha,best.beta);
    int ids[4]={0,1},dim=2;if(best.alpha>0)ids[dim++]=2;if(best.beta>0)ids[dim++]=3;
    double H[4][4]={},cov[4][4]={};bool ok=true;
    for(int i=0;i<dim;++i)for(int j=0;j<dim;++j)H[i][j]=p.H[ids[i]][ids[j]]*c.n;
    for(int j=0;j<dim;++j){double unit[4]={},column[4]={};unit[j]=1;ok=ok&&solve(H,unit,dim,column);for(int i=0;i<dim;++i)cov[ids[i]][ids[j]]=column[i];}
    double root=std::sqrt(best.alpha*best.alpha+best.beta),eq=best.alpha+root;
    if(ok && root>0 && eq>0) {
      double d[4]={0,0,(1+best.alpha/root)/eq,.5/(root*eq)},v=0;
      for(int i=0;i<4;++i)for(int j=0;j<4;++j)v+=d[i]*cov[i][j]*d[j];
      best.seB=std::sqrt(std::max(0.,cov[0][0]))/best.B;best.seEq=std::sqrt(std::max(0.,v));
    }
  }
  return best;
}
}
