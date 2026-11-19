#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <vector>
#include <limits>
#include <stdexcept>
namespace sbnpmle {
using std::vector;
constexpr double pi=3.1415926535897932384626433832795;
inline double normal_pdf(double x){return std::exp(-.5*x*x)/std::sqrt(2*pi);}
inline double normal_cdf(double x){return .5*std::erfc(-x/std::sqrt(2.));}
struct SplineSide {
 vector<double> knots,area,poly; int intervals=0,count=0,active=4;
 explicit SplineSide(vector<double> k,int degree=3):knots(std::move(k)) {
  if(degree==1){
   intervals=(int)knots.size()-1;count=intervals+1;active=2;area.assign(count,0.);
   for(int i=0;i<intervals;++i){double h=knots[i+1]-knots[i];area[i]+=h/2;area[i+1]+=h/2;}
   poly.assign(intervals*16,0.);
   for(int i=0;i<intervals;++i){double h=knots[i+1]-knots[i];poly[i*16]=1/area[i];poly[i*16+1]=-1/(h*area[i]);poly[i*16+5]=1/(h*area[i+1]);}
   return;
  }
  intervals=(int)knots.size()-1; count=intervals+3;
  vector<double> t(intervals+7);
  for(int q=0;q<4;++q){t[q]=0;t[intervals+3+q]=knots.back();}
  for(int i=1;i<intervals;++i)t[3+i]=knots[i];
  area.resize(count);for(int j=0;j<count;++j)area[j]=(t[j+4]-t[j])/4;
  poly.assign(intervals*16,0);
  for(int i=0;i<intervals;++i){
   int mu=i+3;double x0=t[mu],n[4][4]={};n[3][0]=1;
   for(int d=1;d<=3;++d){double nn[4][4]={};
    for(int jj=3-d;jj<=3;++jj){int j=mu-3+jj;double den=t[j+d]-t[j];
     if(den>0)for(int p=0;p<4;++p){nn[jj][p]+=(x0-t[j])/den*n[jj][p];if(p<3)nn[jj][p+1]+=n[jj][p]/den;}
     den=t[j+d+1]-t[j+1];
     if(den>0&&jj<3)for(int p=0;p<4;++p){nn[jj][p]+=(t[j+d+1]-x0)/den*n[jj+1][p];if(p<3)nn[jj][p+1]-=n[jj+1][p]/den;}
    }std::memcpy(n,nn,sizeof n);
   }
   for(int j=0;j<4;++j)for(int p=0;p<4;++p)poly[(i*4+j)*4+p]=n[j][p]/area[i+j];
  }
 }
 void latent(double x,double* pdf,double* cdf) const {
  std::fill(pdf,pdf+count,0);std::fill(cdf,cdf+count,0);
  for(int i=0;i<intervals;++i){double h=knots[i+1]-knots[i],v=std::clamp(x-knots[i],0.,h);
   for(int j=0;j<active;++j){const double*c=&poly[(i*4+j)*4];
    cdf[i+j]+=v*(c[0]+v*(c[1]/2+v*(c[2]/3+v*c[3]/4)));
    if(x>=knots[i]&&(x<knots[i+1]||(i==intervals-1&&x==knots.back())))pdf[i+j]=std::max(0.,c[0]+v*(c[1]+v*(c[2]+v*c[3])));
   }
  }
 }
 void gaussian(double x,double sigma,double* pdf,double* cdf) const {
  if(sigma==0){latent(x,pdf,cdf);return;}
  std::fill(pdf,pdf+count,0);std::fill(cdf,cdf+count,0);
  static constexpr int choose[5][5]={{1,0,0,0,0},{1,1,0,0,0},{1,2,1,0,0},{1,3,3,1,0},{1,4,6,4,1}};
  for(int i=0;i<intervals;++i){double lo=knots[i],hi=knots[i+1],h=hi-lo;
   if(x<lo-9*sigma)continue;
   if(x>hi+9*sigma){for(int j=0;j<active;++j){const double*c=&poly[(i*4+j)*4];cdf[i+j]+=h*(c[0]+h*(c[1]/2+h*(c[2]/3+h*c[3]/4)));}continue;}
   double v0=(lo-x)/sigma,v1=(hi-x)/sigma,p0=normal_pdf(v0),p1=normal_pdf(v1),m[5];
   m[0]=v0>=0 ? normal_cdf(-v0)-normal_cdf(-v1) : normal_cdf(v1)-normal_cdf(v0);
   m[1]=p0-p1;for(int k=2;k<5;++k)m[k]=(k-1)*m[k-2]+std::pow(v0,k-1)*p0-std::pow(v1,k-1)*p1;
   double e[5]={},d=x-lo;
   for(int k=0;k<5;++k)for(int q=0;q<=k;++q)e[k]+=choose[k][q]*std::pow(d,k-q)*std::pow(sigma,q)*m[q];
   double edge=normal_cdf((x-hi)/sigma);
   for(int j=0;j<active;++j){const double*c=&poly[(i*4+j)*4];double f=0,F=0;
    for(int k=0;k<4;++k){f+=c[k]*e[k];F+=c[k]*(std::pow(h,k+1)*edge+e[k+1])/(k+1);}
    pdf[i+j]+=f;cdf[i+j]+=F;
   }
  }
  for(int j=0;j<count;++j){pdf[j]=std::max(0.,pdf[j]);cdf[j]=std::clamp(cdf[j],0.,1.);}
 }
};
inline vector<double> make_knots(double core,double extreme,int intervals){
 vector<double> k{0};double h=core/intervals;
 for(int i=1;i<=intervals;++i)k.push_back(i*h);
 while(k.back()<extreme){h*=2;k.push_back(k.back()+h);if(k.size()>100)throw std::runtime_error("spline support is too large");}
 return k;
}
struct SplineBasis {
 SplineSide left,right;int count;vector<double> area;
 SplineBasis(vector<double> l,vector<double> r,int degree=3):left(std::move(l),degree),right(std::move(r),degree),count(left.count+right.count){area=left.area;area.insert(area.end(),right.area.begin(),right.area.end());}
 void evaluate(double x,double sigma,double* pdf,double* cdf)const{
  left.gaussian(-x,sigma,pdf,cdf);
  for(int j=0;j<left.count;++j)cdf[j]=1-cdf[j];
  right.gaussian(x,sigma,pdf+left.count,cdf+left.count);
  if(sigma==0){if(x>=0)std::fill(pdf,pdf+left.count,0.);else std::fill(pdf+left.count,pdf+count,0.);}
 }
};
// Gauss-Legendre pairs with deterministic subdivision. Reported integration
// errors are numerical estimates; the separate Gaussian tail bound is analytic.
struct FrictionKernel {
 const SplineBasis& basis;double sigma,a,tolerance;int dim;mutable double error=0;mutable bool exhausted=false;mutable bool bin_mode=false;mutable double bin_hi=0;mutable vector<double> component_error;
 FrictionKernel(const SplineBasis& b,double s,double c,double tol):basis(b),sigma(s),a(c),tolerance(tol),dim(b.count+1){}
 void integrand(double v,double x,vector<double>& out)const{
  vector<double> p(dim),c(dim);basis.evaluate(v,sigma,p.data(),c.data());
  p.back()=sigma>0&&std::fabs(v/sigma)<10?normal_pdf(v/sigma)/sigma:0;
  double d=x-v,den=d*d+a*a,k=a/(pi*den),K=.5+std::atan(d/a)/pi;
  if(bin_mode){k=std::atan2(a*(bin_hi-x),a*a+(x-v)*(bin_hi-v))/pi;K=0;}
  out.resize(2*dim);for(int j=0;j<dim;++j){out[j]=p[j]*k;out[dim+j]=p[j]*K;}
 }
 void quadrature(double lo,double hi,double x,bool fine,vector<double>& sum)const{
  static const double x8[]={.183434642495649805,.525532409916328986,.796666477413626740,.960289856497536232};
  static const double w8[]={.362683783378361983,.313706645877887287,.222381034453374471,.101228536290376259};
  static const double x16[]={.0950125098376374402,.281603550779258913,.458016777657227386,.617876244402643748,.755404408355003034,.865631202387831744,.944575023073232576,.989400934991649933};
  static const double w16[]={.189450610455068496,.182603415044923589,.169156519395002538,.149595988816576733,.124628971255533872,.095158511682492785,.062253523938647893,.027152459411754095};
  double mid=(lo+hi)/2,h=(hi-lo)/2;sum.assign(2*dim,0);vector<double> val;
  int n=fine?8:4;for(int q=0;q<n;++q)for(int sign:{-1,1}){integrand(mid+sign*h*(fine?x16[q]:x8[q]),x,val);double w=h*(fine?w16[q]:w8[q]);for(int j=0;j<2*dim;++j)sum[j]+=w*val[j];}
 }
 void integrate(double lo,double hi,double x,double tol,int depth,vector<double>& total)const{
  vector<double> coarse,fine;quadrature(lo,hi,x,false,coarse);quadrature(lo,hi,x,true,fine);
  double err=0,absolute_error=0,scale=std::max(sigma,a);for(int j=0;j<2*dim;++j){double d=std::fabs(coarse[j]-fine[j]);absolute_error=std::max(absolute_error,d);err=std::max(err,bin_mode?d/std::max(1e-24,std::fabs(fine[j])):d*(j<dim?scale:1));}
  if(err<=tol||depth==0){error+=bin_mode?absolute_error:err;for(int j=0;j<2*dim;++j){total[j]+=fine[j];if(bin_mode&&j<dim)component_error[j]+=std::fabs(coarse[j]-fine[j]);}return;}
  double mid=(lo+hi)/2;integrate(lo,mid,x,bin_mode?tol:tol/2,depth-1,total);integrate(mid,hi,x,bin_mode?tol:tol/2,depth-1,total);
 }
 void bin_probability(double lower,double upper,vector<double>& probability)const{
  if(a==0){vector<double> p,l,u;evaluate(lower,p,l);evaluate(upper,p,u);probability.resize(dim);for(int j=0;j<dim;++j)probability[j]=std::max(0.,u[j]-l[j]);return;}
  bin_mode=true;bin_hi=upper;component_error.assign(dim,0.);double lo=-basis.left.knots.back()-9*sigma,hi=basis.right.knots.back()+9*sigma;
  vector<double> cuts{lo,hi,0};for(double v:basis.left.knots)cuts.push_back(-v);for(double v:basis.right.knots)cuts.push_back(v);
  for(double x:{lower,upper})for(double v:{-9*sigma,-sigma,sigma,9*sigma,x-16*a,x-4*a,x-a,x,x+a,x+4*a,x+16*a})if(v>lo&&v<hi)cuts.push_back(v);
  std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());vector<double> sum(2*dim,0.);
  for(size_t i=1;i<cuts.size();++i)integrate(cuts[i-1],cuts[i],lower,tolerance,14,sum);
  probability.assign(sum.begin(),sum.begin()+dim);
  for(int j=0;j<dim;++j)if(component_error[j]>tolerance)exhausted=true;
  if(sigma==0)probability.back()=std::atan2(a*(upper-lower),a*a+upper*lower)/pi;
  bin_mode=false;
 }
 void evaluate(double x,vector<double>& pdf,vector<double>& cdf)const{
  bin_mode=false;double error_before=error;pdf.resize(dim);cdf.resize(dim);
  if(a==0){basis.evaluate(x,sigma,pdf.data(),cdf.data());pdf.back()=normal_pdf(x/sigma)/sigma;cdf.back()=normal_cdf(x/sigma);return;}
  double lo=-basis.left.knots.back()-9*sigma,hi=basis.right.knots.back()+9*sigma;
  vector<double> cuts{lo,hi,0};for(double v:basis.left.knots)cuts.push_back(-v);for(double v:basis.right.knots)cuts.push_back(v);
  for(double v:{-9*sigma,-sigma,sigma,9*sigma,x-16*a,x-4*a,x-a,x,x+a,x+4*a,x+16*a})if(v>lo&&v<hi)cuts.push_back(v);
  std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());vector<double> sum(2*dim,0);
  for(size_t i=1;i<cuts.size();++i)integrate(cuts[i-1],cuts[i],x,tolerance/(cuts.size()-1),14,sum);
  for(int j=0;j<dim;++j){pdf[j]=std::max(0.,sum[j]);cdf[j]=std::clamp(sum[dim+j],0.,1.);}
  if(sigma==0){pdf.back()=a/(pi*(x*x+a*a));cdf.back()=.5+std::atan(x/a)/pi;}
  if(error-error_before>tolerance)exhausted=true;
 }
};
}
