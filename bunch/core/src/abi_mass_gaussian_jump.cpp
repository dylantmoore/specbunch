#include "public_result.hpp"








#include "llcore_prelude.hpp"  
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

namespace cfjv {
#define main cfjv_engine_main_unused
#include "cfj/cfj.cpp"
#undef main
}  

#include "policy_run.hpp"

namespace {

using std::vector;
const double kNaN = std::numeric_limits<double>::quiet_NaN();

struct CfjOut {
  std::string internal;
  double sigma_p = kNaN;
  bool havep = false; int pa = -1, pe = -1; double pB = kNaN, pS = kNaN, pA = kNaN;
  int M2 = 0; double u2 = kNaN;
  vector<double> table;  
  vector<double> u;      
};


static void cfj_core(const vector<double>& dZ, uint64_t seed, CfjOut& O) {
  using namespace cfjv;
  const size_t n = dZ.size();
  Pilot P = ensemble_sigma(dZ);
  O.sigma_p = P.sigma_p;
  if (!(std::isfinite(P.sigma_p) && P.sigma_p > 0)) { O.internal = "no_frequency_band"; return; }
  const double t_start = UGRID_LO / P.sigma_p, Delta = UGRID_DU / P.sigma_p;
  vector<double> Cf, cosine_second_moment;
  cf_recurrence(dZ, t_start, Delta, NGMAX, Cf, cosine_second_moment);
  Ctx cg[3];
  cg[0] = make_ctx(n, P.sigma_p, NG0, Cf, cosine_second_moment);
  vector<double> S, sine_second_moment, CS;
  sin_moments(dZ, t_start, Delta, NGMAX, S, sine_second_moment, CS);
  Mom Mg[3];
  Mg[0] = make_mom(cg[0], S, sine_second_moment, CS);
  vector<double> Yc0, Ys0, Sb, Yc((size_t)BOOT_B * NGPAD, 0.0), Ys((size_t)BOOT_B * NGPAD, 0.0);
  boot_g0<false>(dZ, cg[0].t[0], cg[0].t[1] - cg[0].t[0], seed, Yc0, &Sb);
  for (int b = 0; b < BOOT_B; ++b) std::copy(&Yc0[(size_t)b * MP0], &Yc0[(size_t)b * MP0] + NG0, &Yc[(size_t)b * NGPAD]);
  boot_g0<true>(dZ, cg[0].t[0], cg[0].t[1] - cg[0].t[0], seed, Ys0, nullptr);
  for (int b = 0; b < BOOT_B; ++b) std::copy(&Ys0[(size_t)b * MP0], &Ys0[(size_t)b * MP0] + NG0, &Ys[(size_t)b * NGPAD]);
  Table TJ0 = run_variant(true, cg[0], Mg[0], Yc, Ys, Sb);
  vector<int> rk = s2_rank(TJ0, TJ0.pR0);
  const bool havep = !rk.empty();
  double pB = kNaN, pS = kNaN, pA = kNaN; int pa = -1, pe = -1;
  if (havep) { const VF& o = TJ0.T[rk[0]]; pB = o.r.B; pS = o.r.sigma; pA = o.p1; pa = rk[0] / NG0; pe = rk[0] % NG0; }
  Ctx cmax = make_ctx(n, P.sigma_p, NGMAX, Cf, cosine_second_moment);
  int j2 = NG0 - 1;
  if (havep) {
    for (int j = NG0 - 1; j < NGMAX; ++j) {
      const double t = cmax.t[j], E = std::exp(-pS * pS * t * t / 2.0);
      const double mB = pB * E, mA = pA * E / t;
      const double m = std::sqrt(mB * mB + mA * mA);
      const double Vss = sine_second_moment[j] - S[j] * S[j];
      const double se = std::sqrt((cmax.var[j] + Vss) / (double)n);
      if (m >= se) j2 = j;
    }
  }
  const int M2 = j2 + 1;
  O.havep = havep; O.pa = pa; O.pe = pe; O.pB = pB; O.pS = pS; O.pA = pA; O.M2 = M2; O.u2 = (20.0 + 5.0 * (M2 - 1)) / 100.0;
  const bool g2same = (M2 == NG0);
  Table TG;
  if (g2same) { TG = std::move(TJ0); O.u = cg[0].u; }
  else {
    const int W = std::min(NGPAD, ((std::max(NG1, M2) + 15) / 16) * 16);
    {
      vector<double> Ye;
      boot_cols<false>(dZ, cg[0].t[0], cg[0].t[1] - cg[0].t[0], seed, 64, W, Ye, nullptr);
      for (int b = 0; b < BOOT_B; ++b) std::copy(&Ye[(size_t)b * NGPAD + NG0], &Ye[(size_t)b * NGPAD + W], &Yc[(size_t)b * NGPAD + NG0]);
      Ye.clear();
      boot_cols<true>(dZ, cg[0].t[0], cg[0].t[1] - cg[0].t[0], seed, 64, W, Ye, nullptr);
      for (int b = 0; b < BOOT_B; ++b) std::copy(&Ye[(size_t)b * NGPAD + NG0], &Ye[(size_t)b * NGPAD + W], &Ys[(size_t)b * NGPAD + NG0]);
    }
    cg[2] = make_ctx(n, P.sigma_p, M2, Cf, cosine_second_moment); Mg[2] = make_mom(cg[2], S, sine_second_moment, CS);
    TG = run_variant(true, cg[2], Mg[2], Yc, Ys, Sb);
    O.u = cg[2].u;
  }
  
  vector<int> r2 = s2_rank(TG, TG.pR0);
  if (r2.empty()) O.internal = "no_coherent_window";
  const int NG = TG.M;
  for (int id : r2) {
    const VF& o = TG.T[id]; const BF& r = o.r;
    const bool rok = r.valid && std::isfinite(r.B) && r.B > 0.0 && r.B <= 1.5 && std::isfinite(r.sigma) && r.sigma >= 1e-3 && r.sigma <= 5.0;
    const double row[10] = {(double)(id / NG), (double)(id % NG), r.B, r.sigma, (r.valid && r.sevalid) ? r.seB : kNaN,
                            (r.valid && r.sevalid) ? r.seS : kNaN, TG.pR0[id], rok ? 1.0 : 0.0, kNaN, o.p1};
    O.table.insert(O.table.end(), row, row + 10);
  }
}

int mass_gaussian_jump_impl(const sb_data* data, const sb_opts_mass* opts_in, sb_result** out) {
  if (!data || !out) return sbc::fail(SBC_E_NULL, "data or out is NULL");
  *out = nullptr;
  sb_opts_mass o;
  if (int rc = sbrun::read_opts(opts_in, o)) return rc;
  if (data->w) return sbc::fail(SBC_E_INPUT, "gaussian_jump does not support weights");
  if (data->n == 0) return sbpublic::insufficient(out,0,"gaussian_jump",o.threads);
  if (!data->z) return sbc::fail(SBC_E_INPUT, "sample pointer is NULL");
  vector<double> dZ; dZ.reserve(data->n);
  for (size_t i = 0; i < data->n; ++i) if (std::isfinite(data->z[i])) dZ.push_back(data->z[i]);
  if(dZ.size()<100)return sbpublic::insufficient(out,dZ.size(),"gaussian_jump",o.threads);
  auto* R = new sbc::Result();
  std::unique_ptr<sbc::Result> guard(R);
  CfjOut C;
  vector<double> sums(1 + 2 * 128 + 2 * 66, kNaN);
  int idx = -1;
  {
    sbc::ThreadGuard tg(o.threads);
    R->set("threads_used", tg.count());
    cfj_core(dZ, o.seed, C);
    if (o.distribution_checks) sbdistribution::summarize(dZ, 0.0, sums.data());
    const vector<double>& u = C.u;
    auto uof = [&](int j) { return (j >= 0 && j < (int)u.size()) ? u[j] : kNaN; };
    idx = sbrun::run_policy(R, C.table, 10, C.internal, sums.data(), o, uof, "gaussian_jump", true);
  }
  const size_t nc = C.table.size() / 10;
  R->set("aJ", idx >= 0 ? C.table[(size_t)idx * 10 + 9] : kNaN);
  
  if (!C.internal.empty() || nc == 0) {
    R->set("stage_a.answered", 0); R->sets("stage_a.status", C.internal.empty() ? "no_coherent_window" : C.internal);
    for (const char* k : {"stage_a.B", "stage_a.sigma", "stage_a.aJ", "stage_a.p", "stage_a.a", "stage_a.e"}) R->set(k, kNaN);
  } else {
    const double* q = &C.table[0];
    const bool rok = q[7] != 0.0;
    R->set("stage_a.answered", rok ? 1 : 0); R->sets("stage_a.status", rok ? "ok" : "refused:invalid_fit");
    R->set("stage_a.B", rok ? std::min(std::max(q[2], 0.0), 1.0) : kNaN); R->set("stage_a.sigma", rok ? q[3] : kNaN);
    R->set("stage_a.aJ", rok ? q[9] : kNaN); R->set("stage_a.p", rok ? q[6] : kNaN);
    R->set("stage_a.a", q[0]); R->set("stage_a.e", q[1]);
  }
  R->set("n", (double)dZ.size()); R->set("n_input", (double)data->n); R->set("sigma_pilot", C.sigma_p);
  R->set("pilot.a", C.pa); R->set("pilot.e", C.pe); R->set("pilot.B", C.pB); R->set("pilot.sigma", C.pS); R->set("pilot.aJ", C.pA);
  R->set("M2", C.M2); R->set("u2", C.u2); R->set("seed", (double)o.seed);
  if (true) {
    const char* names[9] = {"cand_a", "cand_e", "cand_B", "cand_sigma", "cand_se_logB", "cand_se_logS", "cand_p", "cand_range_ok", "cand_aJ"};
    const int col[9] = {0, 1, 2, 3, 4, 5, 6, 7, 9};
    for (int f = 0; f < 9; ++f) { vector<double> v(nc); for (size_t i = 0; i < nc; ++i) v[i] = C.table[i * 10 + col[f]]; R->setv(names[f], std::move(v)); }
  }
  if (o.distribution_checks) R->setv("distribution_summaries", sums);
  sbpublic::add_simple(data,o,*R);
  sbpublic::finish(*R,"gaussian_jump");
  R->finish_names();
  *out = reinterpret_cast<sb_result*>(guard.release());
  return SBC_OK;
}

}  

extern "C" {
SB_API int sb_mass_gaussian_jump(const sb_data* data, const sb_opts_mass* opts, sb_result** out) {
  return sbc::guarded([&] { return mass_gaussian_jump_impl(data, opts, out); });
}
}
