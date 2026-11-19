#if defined(__clang__)
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#endif
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <string.h>
#include <math.h>
#include "core/include/specbunch.h"
static void release(SEXP xp) {
  sb_result *r = R_ExternalPtrAddr(xp);
  if (r) { sb_result_free(r); R_ClearExternalPtr(xp); }
}
static SEXP extract(SEXP xp) {
  sb_result *r = R_ExternalPtrAddr(xp);
  int count=0;
  for(int k=0;k<3;k++) {
    const char *s=sb_result_names(r,k);
    if(s && *s) { count++; for(;*s;s++) if(*s=='\n' && s[1]) count++; }
  }
  SEXP raw=PROTECT(Rf_allocVector(VECSXP,count));
  SEXP names=PROTECT(Rf_allocVector(STRSXP,count));
  SEXP kinds=PROTECT(Rf_allocVector(INTSXP,count));
  int i=0;
  for(int k=0;k<3;k++) {
    const char *s=sb_result_names(r,k);
    while(s && *s) {
      const char *end=strchr(s,'\n'); size_t len=end?(size_t)(end-s):strlen(s);
      char *name=(char*)R_alloc(len+1,1); memcpy(name,s,len); name[len]=0;
      SET_STRING_ELT(names,i,Rf_mkChar(name)); INTEGER(kinds)[i]=k;
      if(k==0) { double v; sb_result_scalar(r,name,&v); SET_VECTOR_ELT(raw,i,Rf_ScalarReal(v)); }
      else if(k==1) {
        size_t n=0; const double *v=sb_result_vector(r,name,&n);
        SEXP col=PROTECT(Rf_allocVector(REALSXP,n));
        if(n) memcpy(REAL(col),v,n*sizeof(double)); SET_VECTOR_ELT(raw,i,col); UNPROTECT(1);
      } else SET_VECTOR_ELT(raw,i,Rf_mkString(sb_result_string(r,name)));
      i++; s=end?end+1:NULL;
    }
  }
  Rf_setAttrib(raw,R_NamesSymbol,names); Rf_setAttrib(kinds,R_NamesSymbol,names);
  SEXP out=PROTECT(Rf_allocVector(VECSXP,3));
  SET_VECTOR_ELT(out,0,Rf_ScalarInteger(sb_result_status(r)));
  SET_VECTOR_ELT(out,1,raw); SET_VECTOR_ELT(out,2,kinds);
  release(xp); UNPROTECT(4); return out;
}
static const double *column(SEXP x, R_xlen_t n, int optional) {
  if(optional && x==R_NilValue) return NULL;
  if(TYPEOF(x)!=REALSXP || XLENGTH(x)!=n) Rf_error("Native columns must be double vectors of equal length");
  return REAL(x);
}
static double option(SEXP opts, const char *key, double fallback) {
  SEXP names=Rf_getAttrib(opts,R_NamesSymbol);
  for(R_xlen_t i=0;i<XLENGTH(opts);i++) if(!strcmp(CHAR(STRING_ELT(names,i)),key)) return Rf_asReal(VECTOR_ELT(opts,i));
  return fallback;
}
static SEXP default_threads(void) { return Rf_ScalarInteger(sb_default_threads()); }
static SEXP info(void) {
  SEXP out=PROTECT(Rf_allocVector(VECSXP,3));
  SET_VECTOR_ELT(out,0,Rf_mkString(sb_version()));
  SET_VECTOR_ELT(out,1,Rf_ScalarInteger(sb_abi_version()));
  SET_VECTOR_ELT(out,2,Rf_mkString(sb_build_info()));
  UNPROTECT(1);return out;
}
static SEXP estimate(SEXP entry, SEXP z, SEXP x, SEXP w, SEXP opts) {
  R_xlen_t n=XLENGTH(z);
  const double *zz=column(z,n,0), *ww=column(w,n,1);
  SEXP xp=PROTECT(R_MakeExternalPtr(NULL,R_NilValue,R_NilValue)); R_RegisterCFinalizerEx(xp,release,TRUE);
  sb_result *r=NULL; int rc;
  if(Rf_asInteger(entry)==3) {
    sb_data_xy data={(size_t)n,zz,column(x,n,0),ww}; sb_opts_bmean o; sb_opts_bmean_init(&o);
    o.threads=(int)option(opts,"threads",0); o.seed=(uint64_t)option(opts,"seed",1);
    o.reserved_band=(int)option(opts,"reserved_band",2); o.reserved_visibility=(int)option(opts,"reserved_visibility",1); o.R=(int)option(opts,"R",499);
    rc=sb_bmean(&data,&o,&r);
  } else {
    sb_data data={(size_t)n,zz,ww}; sb_opts_mass o; sb_opts_mass_init(&o);
    o.threads=(int)option(opts,"threads",0); o.seed=(uint64_t)option(opts,"seed",1);
    o.distribution_checks=(int)option(opts,"distribution_checks",1); o.policy=(int)option(opts,"policy",0); o.check_all_windows=(int)option(opts,"check_all_windows",0);
    o.lp_iterations=(int64_t)option(opts,"lp_iterations",5000); o.lp_retry_iterations=(int64_t)option(opts,"lp_retry_iterations",31250);
    int method=Rf_asInteger(entry);
    if(method==0 || method==1 || method==2 || method==4) rc=sb_mass(&data,&o,method==2 || method==4,method==1 || method==4,&r);
    else Rf_error("Unknown native estimator");
  }
  R_SetExternalPtrAddr(xp,r);
  if(rc) { release(xp); Rf_error("Core error %d: %s",rc,sb_last_error()); }
  SEXP out=PROTECT(extract(xp)); UNPROTECT(2); return out;
}
static SEXP npmle_estimate(SEXP z,SEXP w,SEXP bins,SEXP opts) {
 R_xlen_t n=XLENGTH(z),nb=bins==R_NilValue?0:XLENGTH(bins);
 sb_data data={(size_t)n,column(z,n,0),column(w,n,1)};
 const double* edges=nb?column(bins,nb,0):NULL;
 sb_opts_npmle o;sb_opts_npmle_init(&o);
 o.threads=(int)option(opts,"threads",o.threads);
 o.core_intervals=(int)option(opts,"core_intervals",o.core_intervals);
 o.grid_points=(int)option(opts,"grid_points",o.grid_points);
 o.exact_threshold=(int)option(opts,"exact_threshold",o.exact_threshold);
 o.max_iterations=(int)option(opts,"max_iterations",o.max_iterations);
 o.lambda_points=(int)option(opts,"lambda_points",o.lambda_points);
 o.refinement_points=(int)option(opts,"refinement_points",o.refinement_points);
 o.B=option(opts,"B",o.B);
 o.sigma=option(opts,"sigma",o.sigma);
 o.a=option(opts,"a",o.a);
 o.kink=option(opts,"kink",o.kink);
 o.histogram_width=option(opts,"histogram_width",o.histogram_width);
 o.lambda_min=option(opts,"lambda_min",o.lambda_min);
 o.lambda_max=option(opts,"lambda_max",o.lambda_max);
 o.quadrature_tolerance=option(opts,"quadrature_tolerance",o.quadrature_tolerance);
 o.optimality_tolerance=option(opts,"optimality_tolerance",o.optimality_tolerance);
 o.reserved=(int)option(opts,"basis",0);
 SEXP xp=PROTECT(R_MakeExternalPtr(NULL,R_NilValue,R_NilValue));R_RegisterCFinalizerEx(xp,release,TRUE);
 sb_result* r=NULL;int rc=sb_npmle(&data,&o,edges,(size_t)nb,&r);R_SetExternalPtrAddr(xp,r);
 if(rc){release(xp);Rf_error("Core error %d: %s",rc,sb_last_error());}
 SEXP out=PROTECT(extract(xp));UNPROTECT(2);return out;
}
static SEXP post_estimate(SEXP entry, SEXP z, SEXP x, SEXP w, SEXP opts) {
 R_xlen_t n=XLENGTH(z);const double* zz=column(z,n,0);const double* ww=column(w,n,1);
 SEXP xp=PROTECT(R_MakeExternalPtr(NULL,R_NilValue,R_NilValue));R_RegisterCFinalizerEx(xp,release,TRUE);
 sb_result* r=NULL;int rc;
 if(Rf_asInteger(entry)==0){
  sb_opts_nonbuncher_mean o;sb_opts_nonbuncher_mean_init(&o);
  o.threads=(int)option(opts,"threads",o.threads);
  o.degree=(int)option(opts,"degree",o.degree);
  o.B=option(opts,"B",o.B);
  o.sigma=option(opts,"sigma",o.sigma);
  o.a=option(opts,"a",o.a);
  o.kink=option(opts,"kink",o.kink);
  o.mu_b=option(opts,"mu_b",o.mu_b);
  o.se_mu_b=option(opts,"se_mu_b",o.se_mu_b);
  o.bandwidth=option(opts,"bandwidth",o.bandwidth);
  o.bin_width=option(opts,"bin_width",o.bin_width);
  sb_data_xy d={(size_t)n,zz,column(x,n,0),ww};rc=sb_nonbuncher_mean(&d,&o,&r);
 }else{
  sb_opts_eti o;sb_opts_eti_init(&o);
  o.threads=(int)option(opts,"threads",o.threads);
  o.integration_intervals=(int)option(opts,"integration_intervals",o.integration_intervals);
  o.B=option(opts,"B",o.B);
  o.sigma=option(opts,"sigma",o.sigma);
  o.a=option(opts,"a",o.a);
  o.kink=option(opts,"kink",o.kink);
  o.window_lo=option(opts,"window_lo",o.window_lo);
  o.window_hi=option(opts,"window_hi",o.window_hi);
  o.net_below=option(opts,"net_below",o.net_below);
  o.net_above=option(opts,"net_above",o.net_above);
  sb_data d={(size_t)n,zz,ww};rc=Rf_asInteger(entry)==2 ? sb_eti_shoulders(&o,option(opts,"left",NAN),option(opts,"right",NAN),&r) : sb_eti(&d,&o,&r);
 }
 R_SetExternalPtrAddr(xp,r);if(rc){release(xp);Rf_error("Core error %d: %s",rc,sb_last_error());}
 SEXP out=PROTECT(extract(xp));UNPROTECT(2);return out;
}
static const R_CallMethodDef calls[] = {
  {"C_post", (DL_FUNC) &post_estimate, 5},
  {"C_npmle", (DL_FUNC) &npmle_estimate, 4}, {"C_estimate", (DL_FUNC) &estimate, 5}, {"C_info", (DL_FUNC) &info, 0}, {"C_default_threads", (DL_FUNC) &default_threads, 0}, {NULL,NULL,0}
};
void R_init_specbunch(DllInfo *dll) {
  R_registerRoutines(dll,NULL,calls,NULL,NULL);
  R_useDynamicSymbols(dll,FALSE); R_forceSymbols(dll,TRUE);
}
