#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
namespace sbnpmle {
using std::vector;
struct Likelihood {
 int m;size_t n;double effective_n;vector<double> design,atom,counts,penalty,penalty_factor;
 void evaluate(const vector<double>& p,double lambda,double& objective,vector<double>& probabilities,vector<double>* gradient=nullptr,vector<double>* hessian=nullptr)const {
  probabilities.resize(n);double ll=0,pen=0;
  for(size_t i=0;i<n;++i){double g=atom[i];for(int j=0;j<m;++j)g+=design[i*m+j]*p[j];if(!(g>0)){objective=INFINITY;return;}probabilities[i]=g;ll+=counts[i]*std::log(g);}
  for(int j=0;j<m;++j)for(int k=0;k<m;++k)pen+=p[j]*penalty[j*m+k]*p[k];
  objective=(-ll+.5*lambda*pen)/effective_n;
  if(gradient){gradient->assign(m,0);
   #pragma omp parallel for schedule(static)
   for(int j=0;j<m;++j){double g=0;for(size_t i=0;i<n;++i)g-=counts[i]*design[i*m+j]/probabilities[i];for(int k=0;k<m;++k)g+=lambda*penalty[j*m+k]*p[k];(*gradient)[j]=g/effective_n;}
  }
  if(hessian){hessian->assign(m*m,0);
   #pragma omp parallel for schedule(static)
   for(int j=0;j<m;++j)for(int k=0;k<=j;++k){double h=0;for(size_t i=0;i<n;++i)h+=counts[i]*(design[i*m+j]/probabilities[i])*(design[i*m+k]/probabilities[i]);h=(h+lambda*penalty[j*m+k])/effective_n;(*hessian)[j*m+k]=(*hessian)[k*m+j]=h;}
  }
 }
};
inline bool cholesky(vector<double>& h,int m){
 for(int j=0;j<m;++j)for(int k=0;k<=j;++k){double v=h[j*m+k];for(int q=0;q<k;++q)v-=h[j*m+q]*h[k*m+q];if(j==k){if(!(v>0))return false;h[j*m+k]=std::sqrt(v);}else h[j*m+k]=v/h[k*m+k];}return true;
}
inline void solve_cholesky(const vector<double>& l,vector<double>& x){int m=(int)x.size();for(int j=0;j<m;++j){for(int k=0;k<j;++k)x[j]-=l[j*m+k]*x[k];x[j]/=l[j*m+j];}for(int j=m-1;j>=0;--j){for(int k=j+1;k<m;++k)x[j]-=l[k*m+j]*x[k];x[j]/=l[j*m+j];}}
inline vector<double> newton_direction(const vector<double>& h,const vector<double>& g,const vector<double>& p,double mean){
 int m=(int)p.size();vector<int> free;for(int j=0;j<m;++j)if(p[j]>1e-14||g[j]<mean)free.push_back(j);
 vector<double> direction(m,0);
 for(int active_step=0;active_step<m&&!free.empty();++active_step){int nf=(int)free.size();vector<double> scale(nf),matrix(nf*nf),rhs(nf),constraint(nf);
  for(int j=0;j<nf;++j)scale[j]=std::sqrt(std::max(h[free[j]*m+free[j]],1e-30));
  for(int j=0;j<nf;++j){rhs[j]=-g[free[j]]/scale[j];constraint[j]=1/scale[j];for(int k=0;k<nf;++k)matrix[j*nf+k]=h[free[j]*m+free[k]]/(scale[j]*scale[k]);matrix[j*nf+j]+=1e-13;}
  if(!cholesky(matrix,nf))break;solve_cholesky(matrix,rhs);solve_cholesky(matrix,constraint);
  double numerator=0,denominator=0;for(int j=0;j<nf;++j){numerator+=rhs[j]/scale[j];denominator+=constraint[j]/scale[j];}
  double multiplier=numerator/denominator;vector<int> keep;bool blocked=false;
  for(int j=0;j<nf;++j){double d=(rhs[j]-multiplier*constraint[j])/scale[j];if(p[free[j]]<=1e-14&&d<0){blocked=true;direction[free[j]]=0;}else{keep.push_back(free[j]);direction[free[j]]=d;}}
  if(!blocked)return direction;free=keep;
 }
 std::fill(direction.begin(),direction.end(),0.);int best=(int)(std::min_element(g.begin(),g.end())-g.begin());for(int j=0;j<m;++j)direction[j]=-p[j];direction[best]+=1;return direction;
}
struct SplineFit {
 vector<double> p;double lambda=NAN,loglik=NAN,edf=NAN,bic=INFINITY,kkt=INFINITY,edf_ridge=0;bool converged=false;int iterations=0;
};
inline SplineFit optimize(const Likelihood& model,double lambda,const vector<double>& initial,int max_iterations,double tolerance){
 SplineFit fit;fit.lambda=lambda;fit.p=initial;int m=model.m;
 vector<double> prob,gradient,hessian,c(m),y(m),trial(m),temp;
 double objective;
 for(int it=0;it<max_iterations;++it){fit.iterations=it+1;model.evaluate(fit.p,lambda,objective,prob,&gradient,&hessian);
  double mean=0;for(int j=0;j<m;++j)mean+=fit.p[j]*gradient[j];double gap=0,positive=0;
  for(int j=0;j<m;++j){gap=std::max(gap,mean-gradient[j]);positive+=fit.p[j]*std::max(gradient[j]-mean,0.);}
  fit.kkt=std::max(gap,positive);if(fit.kkt<=tolerance){fit.converged=true;break;}
  auto direction=newton_direction(hessian,gradient,fit.p,mean);
  double gd=0,alpha=1;for(int j=0;j<m;++j){gd+=gradient[j]*direction[j];if(direction[j]<0)alpha=std::min(alpha,-fit.p[j]/direction[j]);}
  if(!(gd<0)||!(alpha>0))break;
  double fn=INFINITY;bool accepted=false;
  for(int ls=0;ls<60;++ls){for(int j=0;j<m;++j)trial[j]=std::max(0.,fit.p[j]+alpha*direction[j]);double sum=std::accumulate(trial.begin(),trial.end(),0.);for(double&v:trial)v/=sum;model.evaluate(trial,lambda,fn,temp);
   bool stable=std::isfinite(fn);if(stable)for(size_t i=0;i<prob.size();++i)if(temp[i]<.01*prob[i]){stable=false;break;}
   // Near stationarity, the predicted decrease can be below objective
   // summation roundoff. Accept within that bound; KKT still decides convergence.
   double roundoff=8*std::numeric_limits<double>::epsilon()*std::max(1.,std::fabs(objective));
   if(stable&&fn<=objective+1e-4*alpha*gd+roundoff){accepted=true;break;}alpha*=.5;}
  if(!accepted)break;fit.p.swap(trial);
 }
 model.evaluate(fit.p,lambda,objective,prob,&gradient,&hessian);
 double mean=0;for(int j=0;j<m;++j)mean+=fit.p[j]*gradient[j];double gap=0,pos=0;
 for(int j=0;j<m;++j){gap=std::max(gap,mean-gradient[j]);pos+=fit.p[j]*std::max(gradient[j]-mean,0.);}fit.kkt=std::max(gap,pos);
 // Preserve the iteration target; allow a final certificate margin for objective roundoff.
 fit.converged=fit.kkt<=10*tolerance;
 fit.loglik=0;for(size_t i=0;i<model.n;++i)fit.loglik+=model.counts[i]*std::log(prob[i]);
 // QR of the augmented weighted design evaluates tr(H^-1 I) without
 // subtracting nearly equal matrix entries or forming an inverse.
 vector<double> scale(m,0.),upper(m*m,0.),row(m);
 for(int j=0;j<m;++j){double v=lambda*model.penalty[j*m+j];for(size_t i=0;i<model.n;++i){double a=std::sqrt(model.counts[i])*model.design[i*m+j]/prob[i];v+=a*a;}scale[j]=std::sqrt(v);}
 auto append=[&](vector<double>& v){for(int j=0;j<m;++j){double h=std::hypot(upper[j*m+j],v[j]);if(h==0)continue;double c=upper[j*m+j]/h,s=v[j]/h;upper[j*m+j]=h;for(int k=j+1;k<m;++k){double old=upper[j*m+k];upper[j*m+k]=c*old+s*v[k];v[k]=c*v[k]-s*old;}}};
 for(size_t i=0;i<model.n;++i){for(int j=0;j<m;++j)row[j]=std::sqrt(model.counts[i])*model.design[i*m+j]/prob[i]/scale[j];append(row);}
 for(size_t i=0;i<model.penalty_factor.size()/m;++i){for(int j=0;j<m;++j)row[j]=std::sqrt(lambda)*model.penalty_factor[i*m+j]/scale[j];append(row);}
 vector<double> leverage(model.n,0.);
 #pragma omp parallel for schedule(static)
 for(size_t i=0;i<model.n;++i){vector<double> v(m);for(int j=0;j<m;++j){double a=std::sqrt(model.counts[i])*model.design[i*m+j]/prob[i]/scale[j];for(int k=0;k<j;++k)a-=upper[k*m+j]*v[k];v[j]=a/upper[j*m+j];leverage[i]+=v[j]*v[j];}}
 fit.edf=0;for(double v:leverage)fit.edf+=v;fit.edf_ridge=0;
 fit.bic=-2*fit.loglik+fit.edf*std::log(model.effective_n);
 if(!std::isfinite(fit.bic))fit.converged=false;
 return fit;
}
}
