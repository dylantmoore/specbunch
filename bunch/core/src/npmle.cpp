#include "specbunch.h"
#include "sbc_internal.hpp"
#include "npmle_basis.hpp"
#include "npmle_fit.hpp"
#include <map>
#include <memory>
#include <numeric>
namespace sbnpmle {
struct Row {double lo,hi,count;};
static int run(const sb_data* data,const sb_opts_npmle* supplied,const double* edges,size_t nedges,sb_result** out){
 if(!data||!out)return sbc::fail(SBC_E_NULL,"data or result pointer is NULL");*out=nullptr;
 sb_opts_npmle o;sb_opts_npmle_init(&o);
 if(!supplied)return sbc::fail(SBC_E_INPUT,"NPMLE requires fixed mass and friction estimates");
 if(int rc=sbc::check_opts_header(supplied,8,SB_ABI_VERSION))return rc;
 std::memcpy(&o,supplied,std::min<size_t>(supplied->struct_size,sizeof o));
 if(!(o.B>0&&o.B<1&&o.sigma>=0&&o.a>=0&&o.sigma+o.a>0&&std::isfinite(o.B+o.sigma+o.a+o.kink)))return sbc::fail(SBC_E_INPUT,"mass and friction must be coherent and finite");
 if(o.core_intervals<4||o.core_intervals>64||o.grid_points<25||o.grid_points>10001||o.exact_threshold<0||o.max_iterations<1||o.lambda_points<3||o.lambda_points>101||o.refinement_points<3||o.refinement_points>101||!(o.lambda_min>0&&o.lambda_max>o.lambda_min)||!(o.quadrature_tolerance>0)||!(o.optimality_tolerance>0)||o.histogram_width<0||!std::isfinite(o.histogram_width+o.lambda_min+o.lambda_max+o.quadrature_tolerance+o.optimality_tolerance))return sbc::fail(SBC_E_INPUT,"invalid NPMLE numerical options");
 if(!data->z||data->n==0)return sbc::fail(SBC_E_INPUT,"NPMLE requires a nonempty sample");
 if(nedges&&(nedges<2||!edges))return sbc::fail(SBC_E_INPUT,"bins must contain at least two edges");
 for(size_t i=0;i<nedges;++i)if(!std::isfinite(edges[i])||(i&&edges[i]<=edges[i-1]))return sbc::fail(SBC_E_INPUT,"bin edges must be finite and strictly increasing");
 sbc::ThreadGuard threads(o.threads);
 vector<std::pair<double,double>> sample;sample.reserve(data->n);double total=0,squares=0;
 for(size_t i=0;i<data->n;++i){double z=data->z[i]-o.kink,w=data->w?data->w[i]:1.;if(!std::isfinite(z)||!std::isfinite(w)||w<=0)continue;sample.emplace_back(z,w);total+=w;squares+=w*w;}
 if(sample.size()<2||!(total>0&&std::isfinite(total+squares)))return sbc::fail(SBC_E_INPUT,"NPMLE requires at least two usable observations and finite total weight");
 std::stable_sort(sample.begin(),sample.end(),[](auto&x,auto&y){return x.first<y.first;});
 double neff=total*total/squares,scale=std::max(o.sigma,o.a),width=o.histogram_width>0?o.histogram_width:scale/20;
 double minimum=sample.front().first,maximum=sample.back().first;
 auto quantile=[&](double q){return sample[(size_t)(q*(sample.size()-1))].first;};
 double core_left=std::max(4*scale,-quantile(.005)),core_right=std::max(4*scale,quantile(.995));
 SplineBasis basis(make_knots(core_left,std::max(core_left,-minimum+4*scale),o.core_intervals),make_knots(core_right,std::max(core_right,maximum+4*scale),o.core_intervals));
 int m=basis.count;bool exact=sample.size()<=(size_t)o.exact_threshold;vector<Row> rows;
 if(exact){for(auto&zw:sample)rows.push_back({zw.first,zw.first,zw.second*neff/total});}
 else for(size_t i=0;i<sample.size();){double key=std::floor(sample[i].first/width),w=0;size_t j=i;while(j<sample.size()&&std::floor(sample[j].first/width)==key){w+=sample[j].second;++j;}rows.push_back({key*width,(key+1)*width,w*neff/total});i=j;}
 Likelihood model;model.m=m;model.n=rows.size();model.effective_n=neff;model.design.resize(rows.size()*m);model.atom.resize(rows.size());model.counts.resize(rows.size());model.penalty.assign(m*m,0);
 vector<double> integration_error(rows.size(),0);vector<int> integration_failed(rows.size(),0);
 #pragma omp parallel for schedule(static)
 for(size_t i=0;i<rows.size();++i){FrictionKernel kernel(basis,o.sigma,o.a,o.quadrature_tolerance);vector<double> p,c;if(exact)kernel.evaluate(rows[i].lo,p,c);else kernel.bin_probability(rows[i].lo,rows[i].hi,p);
  for(int j=0;j<m;++j)model.design[i*m+j]=(1-o.B)*p[j];
  model.atom[i]=o.B*p.back();model.counts[i]=rows[i].count;
  integration_error[i]=kernel.error;integration_failed[i]=kernel.exhausted;
 }
 for(int offset:{0,basis.left.count}){int count=offset==0?basis.left.count:basis.right.count;for(int j=0;j<count-2;++j){double d[3]={(1-o.B)/basis.area[offset+j],-2*(1-o.B)/basis.area[offset+j+1],(1-o.B)/basis.area[offset+j+2]};size_t row=model.penalty_factor.size();model.penalty_factor.resize(row+m,0.);for(int u=0;u<3;++u)model.penalty_factor[row+offset+j+u]=d[u];for(int u=0;u<3;++u)for(int v=0;v<3;++v)model.penalty[(offset+j+u)*m+offset+j+v]+=d[u]*d[v];}}
 vector<double> initial(m);double area_sum=std::accumulate(basis.area.begin(),basis.area.end(),0.);for(int j=0;j<m;++j)initial[j]=basis.area[j]/area_sum;
 vector<SplineFit> path;std::map<double,size_t> visited;
 auto fit_lambda=[&](double lambda,const vector<double>& start){auto hit=visited.find(lambda);if(hit!=visited.end())return hit->second;size_t id=path.size();auto fit=optimize(model,lambda,start,o.max_iterations,o.optimality_tolerance);if(!fit.converged){vector<double> restart(m,1./m);auto retry=optimize(model,lambda,restart,2*o.max_iterations,o.optimality_tolerance);if(retry.converged||retry.kkt<fit.kkt)fit=std::move(retry);}path.push_back(std::move(fit));visited[lambda]=id;return id;};
 vector<size_t> coarse;for(int i=0;i<o.lambda_points;++i){double lambda=std::exp(std::log(o.lambda_min)+(std::log(o.lambda_max)-std::log(o.lambda_min))*i/(o.lambda_points-1));size_t id=fit_lambda(lambda,initial);if(coarse.empty()||id!=coarse.back())coarse.push_back(id);initial=path[id].p;}
 auto better=[&](size_t a,size_t b){return path[a].bic<path[b].bic;};
 // Extend the initial logarithmic grid until the minimum is bracketed.
 // Explicit initial endpoints outside the default caps remain valid endpoints.
 const double lower_cap=std::min(1e-6,o.lambda_min),upper_cap=std::max(1e10,o.lambda_max);
 const double log_step=std::max(std::log(10.)/4,(std::log(o.lambda_max)-std::log(o.lambda_min))/(o.lambda_points-1));
 int winner=0;
 for(;;){
  winner=0;for(size_t i=1;i<coarse.size();++i)if(better(coarse[i],coarse[winner]))winner=(int)i;
  if(winner==0 && path[coarse.front()].lambda>lower_cap*(1+1e-12)){
   const auto& endpoint=path[coarse.front()];
   double lambda=std::max(lower_cap,std::exp(std::log(endpoint.lambda)-log_step));
   size_t id=fit_lambda(lambda,endpoint.p);coarse.insert(coarse.begin(),id);
  }else if(winner==(int)coarse.size()-1 && path[coarse.back()].lambda<upper_cap*(1-1e-12)){
   const auto& endpoint=path[coarse.back()];
   double lambda=std::min(upper_cap,std::exp(std::log(endpoint.lambda)+log_step));
   size_t id=fit_lambda(lambda,endpoint.p);coarse.push_back(id);
  }else break;
 }
 size_t best=coarse[winner];bool edge=winner==0||winner==(int)coarse.size()-1;
 if(!edge){double lo=std::log(path[coarse[winner-1]].lambda),hi=std::log(path[coarse[winner+1]].lambda);initial=path[coarse[winner-1]].p;for(int i=0;i<o.refinement_points;++i){size_t id=fit_lambda(std::exp(lo+(hi-lo)*i/(o.refinement_points-1)),initial);initial=path[id].p;if(better(id,best))best=id;}}
 const auto& fit=path[best];edge=fit.lambda<=lower_cap*(1+1e-12)||fit.lambda>=upper_cap*(1-1e-12);auto result=std::make_unique<sbc::Result>();auto& r=*result;
 double quadrature_error=0;bool quadrature_failed=false;for(size_t i=0;i<rows.size();++i){quadrature_error=std::max(quadrature_error,integration_error[i]);quadrature_failed=quadrature_failed||integration_failed[i];}
 r.set("lambda",fit.lambda);r.set("edf",fit.edf);r.set("loglik",fit.loglik);r.set("bic",fit.bic);r.set("lambda_at_grid_edge",edge);r.set("B",o.B);r.set("sigma",o.sigma);r.set("a",o.a);r.set("kink",o.kink);r.set("n",sample.size());r.set("n_input",data->n);r.set("n_eff",neff);r.set("threads_used",threads.count());r.set("histogram_width",width);r.set("exact_observations",exact);r.set("likelihood_rows",rows.size());r.set("basis_count",m);r.set("kkt",fit.kkt);r.set("iterations",fit.iterations);r.set("edf_ridge",fit.edf_ridge);
 vector<double> theta(m);for(int j=0;j<m;++j)theta[j]=(1-o.B)*fit.p[j]/basis.area[j];r.setv("coefficients",theta);r.setv("knots_left",basis.left.knots);r.setv("knots_right",basis.right.knots);
 double left=theta[0],right=theta[basis.left.count];r.set("latent_left",left);r.set("latent_right",right);r.set("latent_jump",right-left);
 auto mixture=[&](const vector<double>& v){double smooth=0;for(int j=0;j<m;++j)smooth+=(1-o.B)*fit.p[j]*v[j];return smooth;};
 vector<double> bins;
 if(nedges){for(size_t i=0;i<nedges;++i)bins.push_back(edges[i]-o.kink);}
 else {
  double h=std::max(scale/4,(core_left+core_right)/4000);bins={minimum,maximum};
  double lo=std::max(minimum,-core_left),hi=std::min(maximum,core_right);
  for(double v=std::ceil(lo/h)*h;v<hi;v+=h)bins.push_back(v);
  for(double v=lo,step=h;v>minimum;step*=1.5){bins.push_back(v);v-=step;}
  for(double v=hi,step=h;v<maximum;step*=1.5){bins.push_back(v);v+=step;}
  std::sort(bins.begin(),bins.end());bins.erase(std::unique(bins.begin(),bins.end()),bins.end());
  if(bins.size()<2)bins={minimum-scale/8,maximum+scale/8};
 }
 size_t nb=bins.size()-1;vector<double> obs(nb,0),fitted(nb),nonbuncher(nb),buncher(nb),latent(nb),binlo(nb),binhi(nb);
 for(auto&zw:sample){auto it=std::upper_bound(bins.begin(),bins.end(),zw.first);if(zw.first==bins.back())obs.back()+=zw.second/total;else if(it!=bins.begin()&&it!=bins.end())obs[(size_t)(it-bins.begin()-1)]+=zw.second/total;}
 vector<double> edge_smooth(bins.size()),edge_atom(bins.size()),edge_err(bins.size());vector<int> edge_fail(bins.size());
 #pragma omp parallel for schedule(static)
 for(size_t i=0;i<bins.size();++i){FrictionKernel kernel(basis,o.sigma,o.a,o.quadrature_tolerance);vector<double> p,c;kernel.evaluate(bins[i],p,c);edge_smooth[i]=mixture(c);edge_atom[i]=o.B*c.back();edge_err[i]=kernel.error;edge_fail[i]=kernel.exhausted;}
 for(size_t i=0;i<nb;++i){binlo[i]=bins[i]+o.kink;binhi[i]=bins[i+1]+o.kink;nonbuncher[i]=std::max(0.,edge_smooth[i+1]-edge_smooth[i]);buncher[i]=std::max(0.,edge_atom[i+1]-edge_atom[i]);fitted[i]=nonbuncher[i]+buncher[i];vector<double> p(m),c(m);basis.evaluate((bins[i]+bins[i+1])/2,0,p.data(),c.data());latent[i]=mixture(p);}
 for(size_t i=0;i<bins.size();++i){quadrature_error=std::max(quadrature_error,edge_err[i]);quadrature_failed=quadrature_failed||edge_fail[i];}
 r.setv("histogram.bin_lo",binlo);r.setv("histogram.bin_hi",binhi);r.setv("histogram.observed_share",obs);r.setv("histogram.fitted_share",fitted);r.setv("histogram.nonbuncher_share",nonbuncher);r.setv("histogram.buncher_share",buncher);r.setv("histogram.latent_counterfactual_density",latent);
 vector<double> grid{minimum,maximum};double gl=std::max(minimum,-core_left),gh=std::min(maximum,core_right),step=(gh-gl)/(o.grid_points-1);
 if(step>0){for(int i=0;i<o.grid_points;++i)grid.push_back(gl+i*step);for(double v=gl,h=step;v>minimum;h*=1.25){grid.push_back(v);v-=h;}for(double v=gh,h=step;v<maximum;h*=1.25){grid.push_back(v);v+=h;}}
 if(minimum<=0&&maximum>=0){grid.push_back(0);if(minimum<0)grid.push_back(std::nextafter(o.kink,-INFINITY)-o.kink);}
 std::sort(grid.begin(),grid.end());grid.erase(std::unique(grid.begin(),grid.end()),grid.end());size_t ng=grid.size();vector<double> lp(ng),lc(ng),op(ng),oc(ng),ec(ng),ge(ng);vector<int> gf(ng);
 #pragma omp parallel for schedule(static)
 for(size_t i=0;i<ng;++i){vector<double> p(m),c(m);basis.evaluate(grid[i],0,p.data(),c.data());lp[i]=mixture(p);lc[i]=mixture(c)+(grid[i]>=0?o.B:0);FrictionKernel kernel(basis,o.sigma,o.a,o.quadrature_tolerance);kernel.evaluate(grid[i],p,c);op[i]=mixture(p)+o.B*p.back();oc[i]=mixture(c)+o.B*c.back();ge[i]=kernel.error;gf[i]=kernel.exhausted;}
 size_t index=0;double running=0;for(size_t i=0;i<ng;++i){while(index<sample.size()&&sample[index].first<=grid[i]){running+=sample[index].second;++index;}ec[i]=running/total;grid[i]+=o.kink;quadrature_error=std::max(quadrature_error,ge[i]);quadrature_failed=quadrature_failed||gf[i];}
 r.setv("curves.grid",grid);r.setv("curves.latent_pdf",lp);r.setv("curves.latent_cdf",lc);r.setv("curves.observed_pdf",op);r.setv("curves.observed_cdf",oc);r.setv("curves.empirical_cdf",ec);
 vector<double> lambdas,lls,edfs,bics,converged,kkts;for(auto&entry:visited){const auto&f=path[entry.second];lambdas.push_back(f.lambda);lls.push_back(f.loglik);edfs.push_back(f.edf);bics.push_back(f.bic);converged.push_back(f.converged);kkts.push_back(f.kkt);}
 r.setv("lambda_path.lambda",lambdas);r.setv("lambda_path.loglik",lls);r.setv("lambda_path.edf",edfs);r.setv("lambda_path.bic",bics);r.setv("lambda_path.converged",converged);r.setv("lambda_path.kkt",kkts);
 r.set("quadrature_error",quadrature_error);r.set("gaussian_tail_bound",o.sigma>0?std::erfc(9/std::sqrt(2.)):0);r.set("gaussian_pdf_error_bound",o.sigma==0?0:o.a>0?std::erfc(9/std::sqrt(2.))/(pi*o.a):normal_pdf(9)/o.sigma);r.set("unconverged_lambdas",std::count(converged.begin(),converged.end(),0.));
 r.status=!fit.converged?1:quadrature_failed?2:0;r.set("status_code",r.status);r.sets("status",!fit.converged?"solver_not_converged":quadrature_failed?"quadrature_not_converged":"converged");r.sets("jump_note","The latent density jump is weakly identified.");r.sets("likelihood",exact?"exact_observations":"full_sample_histogram");r.finish_names();*out=reinterpret_cast<sb_result*>(result.release());return SBC_OK;
}
}
extern "C" {
SB_API void sb_opts_npmle_init(sb_opts_npmle* o){if(!o)return;std::memset(o,0,sizeof *o);o->struct_size=sizeof *o;o->abi_version=SB_ABI_VERSION;o->core_intervals=12;o->grid_points=1001;o->exact_threshold=1000;o->max_iterations=200;o->lambda_points=15;o->refinement_points=9;o->B=NAN;o->sigma=NAN;o->lambda_min=.01;o->lambda_max=10000;o->quadrature_tolerance=1e-9;o->optimality_tolerance=1e-7;}
SB_API int sb_npmle(const sb_data* d,const sb_opts_npmle* o,const double* edges,size_t n,sb_result** out){return sbc::guarded([&]{return sbnpmle::run(d,o,edges,n,out);});}
}
