// Conditional shoulder experiment: known mass/friction, deterministic seeds.
#include "specbunch.h"
#include <omp.h>
#include <random>
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <string>
static double get(sb_result* r,const char* key){double v=NAN;sb_result_scalar(r,key,&v);return v;}
int main(int argc,char** argv){
 int reps=argc>1?std::stoi(argv[1]):100,threads=argc>2?std::stoi(argv[2]):8,cells=argc>4?std::stoi(argv[4]):32;
 std::string file=argc>3?argv[3]:"mc.csv";std::vector<std::string> output(cells*reps);
 auto start=std::chrono::steady_clock::now();omp_set_max_active_levels(1);
 #pragma omp parallel for schedule(dynamic) num_threads(threads)
 for(int task=0;task<cells*reps;++task){
  int cell=task/reps+(argc>5?std::stoi(argv[5]):0),rep=task%reps,law=cell%2,jump=(cell/2)%2,ni=(cell/4)%2,si=cell/8;
  size_t n=ni?1000000:200000;double sigma=std::vector<double>{.02,.03,.05,.08}[si],a=law?.025*sigma:0;
  double e=.3,q=.3,L=.6,R=jump?.9:.6,B=e*q*(L+R)/2,tau=(1-B)/(L+R);
  std::mt19937_64 gen(20261119ULL+cell*100000ULL+rep);std::uniform_real_distribution<double> u(0,1);std::normal_distribution<double> norm(0,1);
  std::vector<double> z(n);for(auto& x:z){double v=u(gen);x=v<B?0:(v<B+L*tau?-1:1)*(-tau*std::log(1-u(gen)));x+=sigma*norm(gen);if(law)x+=a*std::tan(3.14159265358979323846*(u(gen)-.5));}
  sb_data d{n,z.data(),nullptr};sb_opts_eti o;sb_opts_eti_init(&o);o.threads=1;o.B=B;o.sigma=sigma;o.a=a;o.window_lo=2/sigma;o.window_hi=3.5/sigma;o.net_below=std::exp(q);o.net_above=1;
  std::string line;
  for(int basis=-1;basis<2;++basis){sb_result* fit=nullptr;double elapsed=0;auto t=std::chrono::steady_clock::now();int rc=0,status=0;double edge=0,weak=0;
   if(basis<0)rc=sb_eti(&d,&o,&fit);
   else {sb_opts_npmle no;sb_opts_npmle_init(&no);no.threads=1;no.B=B;no.sigma=sigma;no.a=a;no.reserved=basis;no.grid_points=25;no.core_intervals=12;no.histogram_width=sigma/2;no.lambda_points=5;no.refinement_points=3;no.quadrature_tolerance=1e-6;
    sb_result* np=nullptr;rc=sb_npmle(&d,&no,nullptr,0,&np);
    if(!rc){status=sb_result_status(np);edge=get(np,"lambda_at_grid_edge");weak=get(np,"shoulders_weak");rc=sb_eti_shoulders(&o,get(np,"latent_left"),get(np,"latent_right"),&fit);}if(np)sb_result_free(np);
   }
   elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();
   if(rc){std::cerr<<"core error "<<task<<" "<<sb_last_error()<<"\n";status=-rc;}
   for(auto rule:{"saez","continuous"})line+=std::to_string(cell)+","+std::to_string(rep)+","+std::to_string(n)+","+std::to_string(sigma)+","+std::to_string(law)+","+std::to_string(jump)+","+(basis<0?"transform":basis==0?"cubic":"linear")+","+rule+","+std::to_string(fit?get(fit,(std::string(rule)+".e").c_str()):NAN)+","+std::to_string(status)+","+std::to_string(edge)+","+std::to_string(weak)+","+std::to_string(elapsed)+"\n";
   if(fit)sb_result_free(fit);
  }
  output[task]=line;
  #pragma omp critical
  std::cerr<<"completed "<<cell<<":"<<rep<<"\n";
 }
 std::ofstream f(file);f<<"cell,rep,n,sigma,voigt,jump,shoulders,rule,e,status,edge,weak,seconds\n";for(auto& l:output)f<<l;
 std::cerr<<"wall_seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";
}
