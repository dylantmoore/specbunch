#include "private_policy.hpp"


#pragma once
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include "sbc_internal.hpp"
#include "specbunch.h"
#include "window_policy.hpp"
#include "distribution_screen.hpp"
#include "distribution_summaries.hpp"

namespace sbrun {

inline int read_opts(const sb_opts_mass* opts_in, sb_opts_mass& o, bool internal_test=false) {
  sb_opts_mass_init(&o);
  if (!opts_in) return SBC_OK;
  int rc = sbc::check_opts_header(opts_in, 8, SB_ABI_VERSION);
  if (rc) return rc;
  const size_t sz = std::min<size_t>(opts_in->struct_size, sizeof(sb_opts_mass));
  sb_opts_mass tmp; sb_opts_mass_init(&tmp);
  std::memcpy(&tmp, opts_in, sz);
  o = tmp; o.struct_size = sizeof(sb_opts_mass);
  if (o.distribution_checks != 0 && o.distribution_checks != 1) return sbc::fail(SBC_E_INPUT, "distribution_checks must be 0 or 1");
  if ((!internal_test && o.policy!=0 && o.policy!=2) || o.policy < 0 || o.policy > 7) return sbc::fail(SBC_E_INPUT, "window rule must be tiered or first");
  if (o.lp_iterations <= 0) o.lp_iterations = SB_LP_ITER_FIRST > 0 ? SB_LP_ITER_FIRST : (int64_t)1e15;
  if (o.lp_retry_iterations <= 0) o.lp_retry_iterations = SB_LP_ITER_RETRY > 0 ? SB_LP_ITER_RETRY : (int64_t)1e15;
  return SBC_OK;
}


inline std::string flags_str(int f) {
  static const char* N[6] = {"cdf_inconsistent", "moments_inconsistent", "no_valid_distribution", "distribution_check_inconclusive", "residual_cdf_inconsistent", "poor_fit"};
  std::string out;
  for (int k = 0; k < 6; ++k) if (f & (1 << k)) out += (out.empty() ? "" : ",") + std::string(N[k]);
  return out;
}

inline std::string policy_name(const std::string& est, bool distribution_checks, int mode) {
  static const char* names[]={"tiered","strict_test_only","first","fallback_test_only","fallback_inconclusive_test_only","tiered_inconclusive_test_only","distribution_first_test_only","distribution_first_inconclusive_test_only"};
  return est+"."+names[mode]+(distribution_checks?".distribution_checks":".fit_checks");
}




inline int run_policy(sbc::Result* R, const std::vector<double>& table, int cw, const std::string& internal, const double* sums,
                      const sb_opts_mass& o, const std::function<double(int)>& u_of, const std::string& est, bool supports_distribution_checks) {
  using std::vector;
  const double NaN = std::numeric_limits<double>::quiet_NaN();
  const bool distribution_checks = supports_distribution_checks && o.distribution_checks;
  const int mode = o.policy;
  const size_t nc = table.size() / cw;
  vector<sbpol::Cand> cands; cands.reserve(nc);
  for (size_t i = 0; i < nc; ++i) {
    const double* q = &table[i * cw];
    cands.push_back(sbpol::make_cand((int)i + 1, (int)q[0], (int)q[1], q[2], q[3], q[4], q[5], q[6], q[7]));
  }
  std::unique_ptr<sbdistribution::Screen> scr;
  bool distribution_error = false;
  if (distribution_checks) {
    sbdistribution::Summ S = sbdistribution::to_summ(sums);
    if (!sbdistribution::summ_valid(S)) distribution_error = true;  
    else scr.reset(new sbdistribution::Screen(S, o.lp_iterations, o.lp_retry_iterations));
  }
  vector<std::string> distribution_labels(nc);
  vector<int> distribution_evaluated(nc, 0);
  const bool unres_flag_only = mode == SB_POLICY_FALLBACK_NOUNRES || mode == SB_POLICY_TIERED_U || mode == SB_POLICY_DISTRIBUTION_FIRST_INCONCLUSIVE;
  std::function<std::string(int)> cb = nullptr;
  if (scr && mode != SB_POLICY_TOP) cb = [&](int i) {
    if (!distribution_evaluated[i]) { distribution_labels[i] = scr->check(cands[i].B, cands[i].sigma); distribution_evaluated[i] = 1; }
    if (unres_flag_only && distribution_labels[i] == "distribution_check_inconclusive") return std::string();  
    return distribution_labels[i];
  };
  sbpol::PolicyOut P;
  if (false) {
    for (auto& st : P.st) { st.answered = false; st.reason = "error"; }
  } else {
    P = sbpol::policy(cands, internal, cb);  
  }
  const sbpol::Stage& s3 = P.st[3];
  const sbpol::Stage& s1 = P.st[1];
  auto distribution_passes = [&](int i) { return !cb || cb(i).empty(); };
  auto first = [&](bool need_distribution_checks, bool need_fit_test) {  
    for (size_t i = 0; i < nc; ++i) {
      const sbpol::Cand& c = cands[i];
      if (!c.fit_converged || !c.coherent_values) continue;
      if (need_fit_test && !c.fit_test_passed) continue;
      if (need_distribution_checks && !distribution_passes((int)i)) continue;
      return (int)i;
    }
    return -1;
  };
  int idx = -1; std::string reason;
  const bool tiered = mode == SB_POLICY_TIERED || mode == SB_POLICY_TIERED_U || mode == SB_POLICY_DISTRIBUTION_FIRST ||
                      mode == SB_POLICY_DISTRIBUTION_FIRST_INCONCLUSIVE;
  if (mode == SB_POLICY_STRICT) { if (s3.answered) idx = s3.idx; else reason = s3.reason; }
  else if (!s1.answered) reason = s1.reason;                                  
  else if (mode == SB_POLICY_TOP) idx = s1.idx;
  else if (s3.answered) idx = s3.idx;                                         
  else if (!tiered) idx = s1.idx;                                             
  else {
    const bool distribution_first = mode == SB_POLICY_DISTRIBUTION_FIRST || mode == SB_POLICY_DISTRIBUTION_FIRST_INCONCLUSIVE;
    idx = distribution_first ? first(true, false) : first(false, true);
    if (idx < 0) idx = distribution_first ? first(false, true) : first(true, false);
    if (idx < 0) idx = s1.idx;
  }
  int flags = distribution_error?8:0;
  if (idx >= 0) {
    const bool passes_all = mode != SB_POLICY_TOP && cands[idx].fit_test_passed && distribution_passes(idx);
    if (passes_all) { if (scr && distribution_evaluated[idx] && distribution_labels[idx] == "distribution_check_inconclusive") flags |= SB_FLAG_DISTRIBUTION_INCONCLUSIVE; }
    else {
      if (scr) flags |= sbdistribution::check_all(*scr, cands[idx].B, cands[idx].sigma);
      if (!cands[idx].fit_test_passed) flags |= SB_FLAG_POOR_FIT;
    }
  }
  const double Bout = idx >= 0 ? cands[idx].B : NaN;  
  const int distribution_flags = flags & (SB_FLAG_CDF_INCONSISTENT | SB_FLAG_MOMENTS_INCONSISTENT | SB_FLAG_NO_VALID_DISTRIBUTION | SB_FLAG_DISTRIBUTION_INCONCLUSIVE | SB_FLAG_RESIDUAL_CDF_INCONSISTENT);
  const int tier = idx < 0 ? 0 : (!distribution_flags && !(flags & SB_FLAG_POOR_FIT)) ? 1 : (distribution_flags && !(flags & SB_FLAG_POOR_FIT)) ? 2 : (!distribution_flags) ? 3 : 4;
  R->set("tier", tier);
  R->sets("policy", policy_name(est, distribution_checks, mode));
  R->set("policy_mode", mode);
  R->sets("fit_failure", internal);
  R->set("n_cands", (double)nc); R->set("n_screened", P.n_screened);
  R->sets("higher", idx >= 0 ? sbpol::rle(P.labels, (size_t)idx) : (mode == SB_POLICY_STRICT ? s3.higher : s1.higher));
  if (idx >= 0) {
    const sbpol::Cand& c = cands[idx];
    R->status = SBC_STATUS_ESTIMATE;
    R->set("B", Bout); R->set("sigma", c.sigma); R->set("se_logB", c.seB); R->set("se_logS", c.seS); R->set("p", c.p);
    R->set("a", c.a); R->set("e", c.e); R->set("u_lo", u_of(c.a)); R->set("u_hi", u_of(c.e)); R->set("rank", c.rank);
    R->sets("win", c.win); R->sets("reason", "");
  } else {
    R->status = reason == "no_coherent_window" ? (P.st[0].answered ? SB_REFUSED_INCOHERENT : SB_REFUSED_NO_FIT) : reason == "no_window_passes:poor_fit" ? SB_REFUSED_POOR_FIT
              : reason.rfind("no_window_passes:check_distribution", 0) == 0 ? SB_REFUSED_DISTRIBUTION : SB_REFUSED_NO_FIT;
    for (const char* k : {"B", "sigma", "se_logB", "se_logS", "p", "a", "e", "u_lo", "u_hi", "rank"}) R->set(k, NaN);
    R->sets("win", ""); R->sets("reason", reason);
  }
  R->set("flags", flags); R->sets("flags_str", flags_str(flags)); R->set("fallback", flags != 0 ? 1.0 : 0.0);
  R->set("status", R->status);
  R->sets("se_note", "se_logB / se_logS are model-based standard errors of the selected fit; they are NOT a confidence interval");
  const char* nm[4] = {"stage.fit_converged", "stage.coherent_values", "stage.distribution_checks", "stage.fit_test"};
  for (int k = 0; k < 4; ++k) {
    R->set(std::string(nm[k]) + ".answered", P.st[k].answered ? 1.0 : 0.0);
    R->set(std::string(nm[k]) + ".rank", P.st[k].answered ? P.st[k].rank : NaN);
    R->sets(std::string(nm[k]) + ".reason", P.st[k].answered ? "" : P.st[k].reason);
  }
  std::string tr;
  for (size_t i = 0; i < nc; ++i)
    if (distribution_evaluated[i]) tr += (tr.empty() ? "" : ";") + std::to_string(i + 1) + ":" + (distribution_labels[i].empty() ? "pass" : distribution_labels[i]);
  R->sets("distribution_trace", tr);
  R->set("distribution_checks", distribution_checks ? 1 : 0); R->set("lp_iterations", (double)o.lp_iterations); R->set("lp_retry_iterations", (double)o.lp_retry_iterations);
  R->set("distribution_lp_iterations", scr ? (double)scr->lp_iters : 0.0);
  R->set("distribution_checks", scr ? scr->n_checks : 0); R->set("distribution_retries", scr ? scr->n_retry : 0);
  R->set("distribution_unresolved", scr ? scr->n_unresolved : 0);
  R->set("distribution_max_lp_iters_resolved", scr ? (double)scr->engine.max_lp_iters_resolved : 0.0);
  R->set("distribution_max_call_iters_first", scr ? (double)scr->max_call_first : 0.0);
  R->set("distribution_max_call_iters_retry", scr ? (double)scr->max_call_retry : 0.0);
  const char* check_names[]={"cdf_inconsistent","moments_inconsistent","no_valid_distribution","distribution_check_inconclusive","residual_cdf_inconsistent"};
  vector<std::array<int,5>> outcomes(nc,{-1,-1,-1,-1,-1});
  if(scr) for(size_t i=0;i<nc;++i) {
    auto it=scr->outcomes.find({cands[i].B,cands[i].sigma});
    if(distribution_evaluated[i] && it!=scr->outcomes.end()) outcomes[i]=it->second;
  }
  if(scr && idx>=0) {
    for(int j=0;j<5;++j) outcomes[idx][j]=(flags&(1<<j))?0:1;
    if(flags&2) outcomes[idx][2]=outcomes[idx][3]=-1;
  }
  if(scr && o.check_all_windows==2) for(size_t i=0;i<nc;++i) if(cands[i].fit_converged&&cands[i].coherent_values && (int)i!=idx) {
    int f=sbdistribution::check_all(*scr,cands[i].B,cands[i].sigma);
    for(int j=0;j<5;++j) outcomes[i][j]=(f&(1<<j))?0:1;
    if(f&2) outcomes[i][2]=outcomes[i][3]=-1;
  }
  for(int j=0;j<5;++j) {vector<double> v(nc);for(size_t i=0;i<nc;++i)v[i]=outcomes[i][j];R->setv(std::string("windows.")+check_names[j],std::move(v));}
  for(const std::string name:{"rank","window_lo_scaled","window_hi_scaled","n_points","B","sigma","se_log_B","se_log_sigma","fit_p","fit_converged","coherent_values","coherent","poor_fit","selected"}) {
    vector<double> v(nc);
    for(size_t i=0;i<nc;++i) {const auto &c=cands[i];
      if(name=="rank")v[i]=i+1;else if(name=="window_lo_scaled")v[i]=u_of(c.a);else if(name=="window_hi_scaled")v[i]=u_of(c.e);
      else if(name=="n_points")v[i]=c.e-c.a+1;else if(name=="B")v[i]=c.B;else if(name=="sigma")v[i]=c.sigma;
      else if(name=="se_log_B")v[i]=c.seB;else if(name=="se_log_sigma")v[i]=c.seS;else if(name=="fit_p")v[i]=c.p;
      else if(name=="fit_converged")v[i]=c.fit_converged?1:0;else if(name=="poor_fit")v[i]=std::isfinite(c.p)?(c.fit_test_passed?1:0):-1;
      else if(name=="selected")v[i]=(int)i==idx;else v[i]=c.fit_converged&&c.coherent_values;
    }R->setv("windows."+name,std::move(v));
  }
  return idx;
}

}  
