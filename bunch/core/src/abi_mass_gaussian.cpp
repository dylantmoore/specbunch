#include "private_policy.hpp"






#include <functional>
#include <memory>
#include "llcore_prelude.hpp"
#include "llcore.hpp"
#include "p1lad.hpp"
#include "p1w.hpp"
#include "puf/p1famv.hpp"
#include "sbc_internal.hpp"
#include "specbunch.h"
#include "public_result.hpp"
#include "window_policy.hpp"
#include "distribution_screen.hpp"
#include "distribution_summaries.hpp"
#include "policy_run.hpp"

namespace {


enum St { ST_OK = 0, ST_PILOT, ST_PILOTB, ST_UH, ST_NOBANDS, ST_NOVALID, ST_INVALID, ST_FEWPOS };
const char* const STNAME[8] = {"ok", "no_frequency_band", "no_frequency_band", "no_coherent_window", "no_coherent_window",
                               "no_coherent_window", "no_coherent_window", "insufficient_data"};
const char* const EST[6] = {"quadratic_reference", "simple_ols", "simple_lad", "ols_ranked", "lad_ranked", "gaussian.rank1"};
// Vendor policy ordinal 1 selects the fixed reference band.
constexpr Policy reference_band_policy = static_cast<Policy>(1);
constexpr int NEST = 6, CW = 9;

struct Est { int st = ST_NOVALID; double B = NaN, sigma = NaN; int a = -1, e = -1, m = 0, ndrop = 0; double p = NaN; };
static void fin(Est& o, const BF& r) {
  Pick pk; finalize(pk, r);
  o.st = pk.status == "ok" ? ST_OK : ST_INVALID; o.B = pk.B; o.sigma = pk.sigma;
}

struct Spec {
  double n_eff = NaN, W = NaN, sigma_p = NaN, sig_filt = NaN, sig_iter = NaN, B_p = NaN, uH = NaN;
  double nv[3] = {0, 0, 0};
  Est est[NEST];
  int fit_failure = ST_OK;
  vector<double> cands;  
  bool have_ctx = false;
  vector<double> u;      
};


static void spectral_core(const vector<double>& dZ, const vector<double>* w, uint64_t seed, Spec& O, bool want_cands) {
  const size_t n = dZ.size();
  p1w::WS ws;
  if (w) { ws = p1w::wsums(*w); O.n_eff = ws.neff; O.W = ws.W; }
  else { O.n_eff = (double)n; O.W = (double)n; }
  for (auto& e : O.est) e = Est();
  O.cands.clear(); O.fit_failure = ST_OK;
  Pilot P = w ? p1w::ensemble_sigma_w(dZ, *w, ws) : ensemble_sigma(dZ);
  O.sigma_p = P.sigma_p; O.sig_filt = P.sig_filt; O.sig_iter = P.sig_iter;
  if (!(std::isfinite(P.sigma_p) && P.sigma_p > 0)) {
    for (auto& e : O.est) e.st = ST_PILOT;
    O.fit_failure = ST_PILOT;
    return;
  }
  Ctx c = w ? p1w::ugrid_ctx_w(dZ, *w, ws, P.sigma_p) : ugrid_ctx(dZ, P.sigma_p);
  O.have_ctx = true; O.u = c.u;
  const double Bp = pilot_B(c);
  vector<double> ssm = running_median5(c.s);
  const int uHs = (std::isfinite(Bp) && Bp > 0) ? choose_uH_idx(reference_band_policy, c, Bp, ssm) : -1;
  O.B_p = Bp; O.uH = uHs;
  vector<double> xs(NG);
  for (int j = 0; j < NG; ++j) xs[j] = -c.t[j] * c.t[j] / 2.0;
  
  vector<BF> T[3];
  vector<double> MU[3];
  for (int G = 0; G < 3; ++G) { T[G].assign(NG * NG, BF()); MU[G].assign((size_t)NG * NG * NG, 0.0); }
#pragma omp parallel for schedule(dynamic, 1)
  for (int a = 0; a < NG; ++a)
    for (int e = a + MINPTS - 1; e < NG; ++e) {
      const int id = bidx(a, e);
      for (int G = 0; G < 3; ++G) {
        BF r;
        bool v = (G == 0) ? fit_band<true>(0, c, xs.data(), c.C.data(), a, e - a + 1, P.sigma_p, r)
               : (G == 1) ? fit_lad(&c.C[a], &xs[a], e - a + 1, NaN, NaN, r)
                          : fit_band<true>(2, c, xs.data(), c.C.data(), a, e - a + 1, P.sigma_p, r);
        if (v) {
          r.D = sup_std_resid(&c.C[a], &c.t[a], &c.s[a], e - a + 1, r.B, r.sigma);
          double* mu = &MU[G][(size_t)id * NG];
          for (int j = a; j <= e; ++j) mu[j] = mu_k0(r.B, r.sigma, c.t[j]);
        }
        r.valid = v;
        T[G][id] = r;
      }
    }
  
  vector<double> Y, S;
  if (w) p1w::boot_sums_w(dZ, *w, c.t[0], c.t[1] - c.t[0], seed, Y, S);
  else boot_sums(dZ, c.t[0], c.t[1] - c.t[0], seed, Y, S);
  const double inv = w ? 1.0 / ws.W : 1.0 / (double)n;
  vector<int> BL[3];
  for (int G = 0; G < 3; ++G)
    for (int a = 0; a < NG; ++a)
      for (int e = a + MINPTS - 1; e < NG; ++e)
        if (T[G][bidx(a, e)].valid && std::isfinite(T[G][bidx(a, e)].D)) BL[G].push_back(bidx(a, e));
  for (int G = 0; G < 3; ++G) O.nv[G] = (double)BL[G].size();
  vector<std::pair<int, int>> tasks;
  for (int G = 0; G < 3; ++G) for (int id : BL[G]) tasks.push_back({G, id});
  const int nt = (int)tasks.size();
  vector<int> counts(nt, 0);
#pragma omp parallel
  {
    vector<int> cl(nt, 0);
    double Cs[NG];
#pragma omp for schedule(dynamic, 4)
    for (int b = 0; b < BOOT_B; ++b) {
      const double* yb = &Y[(size_t)b * MP];
      const double Sb = S[b];
      for (int q = 0; q < nt; ++q) {
        const int G = tasks[q].first, id = tasks[q].second, a = id / NG, e = id % NG, m = e - a + 1;
        const BF* r0 = &T[G][id];
        const double* mub = &MU[G][(size_t)id * NG];
        for (int j = a; j <= e; ++j) Cs[j] = mub[j] + (yb[j] - c.C[j] * Sb) * inv;
        BF r; bool ok;
        if (G == 0) ok = fit_band<false>(0, c, xs.data(), Cs, a, m, r0->sigma, r);
        else if (G == 1) ok = fit_lad(&Cs[a], &xs[a], m, r0->logB, r0->sigma * r0->sigma, r);
        else ok = fit_band<false>(2, c, xs.data(), Cs, a, m, r0->sigma, r);
        if (!ok) continue;
        double Ds = sup_std_resid(&Cs[a], &c.t[a], &c.s[a], m, r.B, r.sigma);
        if (Ds >= r0->D) cl[q]++;
      }
    }
#pragma omp critical
    for (int q = 0; q < nt; ++q) counts[q] += cl[q];  
  }
  vector<double> PV[3];
  for (int G = 0; G < 3; ++G) PV[G].assign(NG * NG, NaN);
  for (int q = 0; q < nt; ++q) PV[tasks[q].first][tasks[q].second] = (1.0 + counts[q]) / (BOOT_B + 1.0);
  
  int gq_st = ST_OK; vector<int> gq_lo;
  if (!(std::isfinite(Bp) && Bp > 0)) gq_st = ST_PILOTB;
  else if (!(uHs > 0 && uHs <= NG)) gq_st = ST_UH;
  else { gq_lo = candidate_lo(uHs, c.u); if (gq_lo.empty()) gq_st = ST_NOBANDS; }
  if (gq_st != ST_OK) { O.est[0].st = O.est[1].st = O.est[2].st = gq_st; }
  else {
    const int eG = uHs - 1;
    vector<int> hi, rest;
    for (int j : gq_lo) {
      const int bid = bidx(j - 1, eG); const BF& r = T[2][bid];
      if (!(r.valid && std::isfinite(r.D))) continue;
      (PV[2][bid] >= 0.10 ? hi : rest).push_back(bid);
    }
    std::stable_sort(rest.begin(), rest.end(), [&](int x, int y) { return PV[2][x] > PV[2][y]; });
    vector<int> rk(hi); rk.insert(rk.end(), rest.begin(), rest.end());
    if (rk.empty()) { O.est[0].st = O.est[1].st = O.est[2].st = ST_NOVALID; }
    else {
      const int bid = rk[0], a = bid / NG, e = eG, m = e - a + 1;
      Est& g = O.est[0]; g.a = a; g.e = e; g.m = m; g.p = PV[2][bid];
      fin(g, T[2][bid]);
      double xx[NG], CC[NG], VV[NG]; int mm = 0;
      for (int j = a; j <= e; ++j) if (c.C[j] > 0.0) { xx[mm] = xs[j]; CC[mm] = c.C[j]; VV[mm] = c.var[j]; ++mm; }
      for (int G = 0; G < 2; ++G) {
        Est& o = O.est[1 + G];
        o.a = a; o.e = e; o.m = m; o.ndrop = m - mm;
        if (mm < 3) { o.st = ST_FEWPOS; continue; }
        BF r;
        bool v = (G == 0) ? fit_log<false, false>(CC, xx, VV, mm, (double)c.n, r) : fit_lad(CC, xx, mm, NaN, NaN, r);
        if (!v) { o.st = ST_INVALID; continue; }
        fin(o, r);
      }
    }
  }
  
  auto ranked_window_order = [&](int G) {
    vector<int> hi, lo;
    for (int id : BL[G]) (PV[G][id] >= 0.10 ? hi : lo).push_back(id);
    auto pts = [](int id) { return id % NG - id / NG + 1; };
    std::stable_sort(hi.begin(), hi.end(), [&](int x, int y) {
      if (pts(x) != pts(y)) return pts(x) > pts(y);
      return x / NG < y / NG; });
    std::stable_sort(lo.begin(), lo.end(), [&](int x, int y) {
      if (PV[G][x] != PV[G][y]) return PV[G][x] > PV[G][y];
      if (pts(x) != pts(y)) return pts(x) > pts(y);
      return x / NG < y / NG; });
    vector<int> rk(hi); rk.insert(rk.end(), lo.begin(), lo.end());
    return rk;
  };
  for (int G = 0; G < 3; ++G) {
    Est& o = O.est[G == 2 ? 5 : 3 + G];
    vector<int> rk = ranked_window_order(G);
    if (rk.empty()) { o.st = ST_NOVALID; if (G == 2) O.fit_failure = ST_NOVALID; continue; }
    const int id = rk[0];
    o.a = id / NG; o.e = id % NG; o.m = o.e - o.a + 1; o.p = PV[G][id]; o.ndrop = 0;
    fin(o, T[G][id]);
    if (G == 2) {  
      (void)want_cands;
      O.cands.reserve(rk.size() * CW);
      for (int q : rk) {
        const BF& r = T[2][q];
        const bool rok = std::isfinite(r.B) && r.B > 0.0 && r.B <= 1.5 && std::isfinite(r.sigma) && r.sigma >= 1e-3 && r.sigma <= 5.0;
        const double row[CW] = {(double)(q / NG), (double)(q % NG), r.B, r.sigma, r.sevalid ? r.seB : NaN, r.sevalid ? r.seS : NaN,
                                PV[2][q], rok ? 1.0 : 0.0, NaN};
        O.cands.insert(O.cands.end(), row, row + CW);
      }
    }
  }
}


static double u_of(int j) {
  return (j == UGRID_N - 1) ? UGRID_HI : UGRID_LO + j * ((UGRID_HI - UGRID_LO) / (UGRID_N - 1));
}

int mass_gaussian_impl(const sb_data* data, const sb_opts_mass* opts_in, sb_result** out) {
  if (!data || !out) return sbc::fail(SBC_E_NULL, "data or out is NULL");
  *out = nullptr;
  sb_opts_mass o;
  if (int rc = sbrun::read_opts(opts_in, o)) return rc;
  if(data->n==0)return sbpublic::insufficient(out,0,"gaussian",o.threads);
  if(!data->z)return sbc::fail(SBC_E_INPUT,"z is NULL");
  
  const double* z = data->z; const double* w = data->w; const size_t n = data->n;
  vector<double> dZ, wv; dZ.reserve(n); if (w) wv.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    if (!std::isfinite(z[i])) continue;
    if (w && !(std::isfinite(w[i]) && w[i] > 0)) continue;
    dZ.push_back(z[i]); if (w) wv.push_back(w[i]);
  }
  if(dZ.size()<100)return sbpublic::insufficient(out,dZ.size(),"gaussian",o.threads);

  auto* R = new sbc::Result();
  std::unique_ptr<sbc::Result> guard(R);
  Spec O;
  vector<double> sums(1 + 2 * 128 + 2 * 66, NaN);
  {
    sbc::ThreadGuard tg(o.threads);
    R->set("threads_used", tg.count());
    spectral_core(dZ, w ? &wv : nullptr, o.seed, O, o.check_all_windows != 0);
    if (o.distribution_checks) { if (w) sbdistribution::summarize_w(dZ, wv, sums.data()); else sbdistribution::summarize(dZ, 0.0, sums.data()); }
    const std::string internal = O.fit_failure == ST_OK ? std::string() : std::string(STNAME[O.fit_failure]);
    sbrun::run_policy(R, O.cands, CW, internal, sums.data(), o, u_of, "gaussian", true);
  }
  R->set("n", (double)dZ.size()); R->set("n_input", (double)n); R->set("n_eff", O.n_eff); R->set("W", O.W);
  R->set("weighted", w ? 1.0 : 0.0);
  R->set("sigma_pilot", O.sigma_p); R->set("sig_filt", O.sig_filt); R->set("sig_iter", O.sig_iter);
  R->set("B_pilot", O.B_p); R->set("reference_band_hi", O.uH);
  R->set("n_valid_ols", O.nv[0]); R->set("n_valid_lad", O.nv[1]); R->set("n_valid_f3", O.nv[2]);
  R->set("seed", (double)o.seed);
  for (int k = 0; k < NEST; ++k) {
    const Est& e = O.est[k]; const std::string p = EST[k];
    R->set(p + ".status", e.st); R->sets(p + ".status_name", STNAME[e.st]);
    R->set(p + ".B", e.B); R->set(p + ".sigma", e.sigma);
    R->set(p + ".uL", (O.have_ctx && e.a >= 0) ? O.u[e.a] : NaN); R->set(p + ".uH", (O.have_ctx && e.e >= 0) ? O.u[e.e] : NaN);
    R->set(p + ".m", e.m); R->set(p + ".ndrop", e.ndrop); R->set(p + ".p", e.p);
  }
  const size_t nc = O.cands.size() / CW;
  if (true) {
    const char* names[8] = {"cand_a", "cand_e", "cand_B", "cand_sigma", "cand_se_logB", "cand_se_logS", "cand_p", "cand_range_ok"};
    for (int f = 0; f < 8; ++f) {
      vector<double> v(nc);
      for (size_t i = 0; i < nc; ++i) v[i] = O.cands[i * CW + f];
      R->setv(names[f], std::move(v));
    }
  }
  if (o.distribution_checks) R->setv("distribution_summaries", sums);
  sbpublic::finish(*R,"gaussian");
  R->finish_names();
  *out = reinterpret_cast<sb_result*>(guard.release());
  return SBC_OK;
}

int policy_impl(const double* table8, size_t ncand, const char* internal, const double* sums, const sb_opts_mass* opts_in,
                sb_result** out) {
  if (!out || (ncand > 0 && !table8)) return sbc::fail(SBC_E_NULL, "NULL argument");
  *out = nullptr;
  sb_opts_mass o;
  if (int rc = sbrun::read_opts(opts_in, o, true)) return rc;
  if (o.distribution_checks && !sums) return sbc::fail(SBC_E_NULL, "summaries are NULL with distribution_checks = 1");
  vector<double> table(ncand * CW, NaN);
  for (size_t i = 0; i < ncand; ++i) for (int f = 0; f < 8; ++f) table[i * CW + f] = table8[i * 8 + f];
  auto* R = new sbc::Result();
  std::unique_ptr<sbc::Result> guard(R);
  {
    sbc::ThreadGuard tg(o.threads);
    R->set("threads_used", tg.count());
    sbrun::run_policy(R, table, CW, internal ? std::string(internal) : std::string(), sums, o, u_of, "gaussian", true);
  }
  R->finish_names();
  *out = reinterpret_cast<sb_result*>(guard.release());
  return SBC_OK;
}



static int filter_rows(const sb_data* data, vector<double>& dZ, vector<double>& wv) {
  const double* z = data->z; const double* w = data->w; const size_t n = data->n;
  dZ.reserve(n); if (w) wv.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    if (!std::isfinite(z[i])) continue;
    if (w && !(std::isfinite(w[i]) && w[i] > 0)) continue;
    dZ.push_back(z[i]); if (w) wv.push_back(w[i]);
  }

  return SBC_OK;
}
static double u_ext(int j) { return 0.20 + 0.05 * j; }  
static const char* FAMV_INTERNAL[3] = {"", "refused:pilot", "refused:no_valid_band"};

int mass_voigt_impl(const sb_data* data, const sb_opts_mass* opts_in, sb_result** out, int fam) {
  if (!data || !out) return sbc::fail(SBC_E_NULL, "data or out is NULL");
  *out = nullptr;
  sb_opts_mass o;
  if (int rc = sbrun::read_opts(opts_in, o)) return rc;
  if(data->n==0)return sbpublic::insufficient(out,0,"voigt",o.threads);
  if(!data->z)return sbc::fail(SBC_E_INPUT,"z is NULL");
  vector<double> dZ, wv;
  if (int rc = filter_rows(data, dZ, wv)) return rc;
  if(dZ.size()<100)return sbpublic::insufficient(out,dZ.size(),"voigt",o.threads);
  const bool weighted = data->w != nullptr;
  auto* R = new sbc::Result();
  std::unique_ptr<sbc::Result> guard(R);
  p1v::Res V;
  sb_opts_mass onc = o; onc.distribution_checks = 0;  
                                     
  int idx = -1;
  {
    sbc::ThreadGuard tg(o.threads);
    R->set("threads_used", tg.count());
    p1v::famv_core(dZ, weighted ? &wv : nullptr, o.seed, fam, V);
    const std::string intl = (V.internal >= 0 && V.internal <= 2) ? FAMV_INTERNAL[V.internal] : "no_coherent_window";
    idx = sbrun::run_policy(R, V.cands, 9, intl, nullptr, onc, u_ext, fam == 2 ? "voigt" : "famlev.F" + std::to_string(fam), false);
  }
  const size_t nc = V.cands.size() / 9;
  
  R->set("sigma_eq", R->scalars["sigma"]);
  if (idx >= 0) {
    const double a = V.x[(size_t)idx * 3], s2 = V.x[(size_t)idx * 3 + 1];
    R->set("a_voigt", a); R->set("s2_voigt", s2); R->set("sigma", s2 >= 0 ? std::sqrt(s2) : NaN);
    R->set("kind", V.x[(size_t)idx * 3 + 2]); R->set("se_log_sigma_eq", R->scalars["se_logS"]);
  } else {
    for (const char* k : {"a_voigt", "s2_voigt", "kind", "se_log_sigma_eq"}) R->set(k, NaN);
  }
  
  if (nc > 0 && V.internal == 0) {
    const double* q = &V.cands[0];
    R->set("rank1.B", q[2]); R->set("rank1.sigma_eq", q[3]); R->set("rank1.a_voigt", V.x[0]);
    R->set("rank1.sigma", V.x[1] >= 0 ? std::sqrt(V.x[1]) : NaN); R->set("rank1.p", q[6]); R->set("rank1.range_ok", q[7]);
    R->set("rank1.u_lo", u_ext((int)q[0])); R->set("rank1.u_hi", u_ext((int)q[1]));
  } else {
    for (const char* k : {"rank1.B", "rank1.sigma_eq", "rank1.a_voigt", "rank1.sigma", "rank1.p", "rank1.range_ok", "rank1.u_lo", "rank1.u_hi"}) R->set(k, NaN);
  }
  R->set("n", (double)dZ.size()); R->set("n_input", (double)data->n); R->set("weighted", weighted ? 1.0 : 0.0);
  R->set("sigma_pilot", V.sigma_p); R->set("NGE", V.NGE); R->set("u_hi_ext", V.u_hi_ext); R->set("n_bands", V.n_bands);
  R->set("kkt_restarts", V.kkt_restarts); R->set("seed", (double)o.seed); R->set("family", fam);
  if (true) {
    const char* names[8] = {"cand_a", "cand_e", "cand_B", "cand_sigma_eq", "cand_se_logB", "cand_se_log_sigma_eq", "cand_p", "cand_range_ok"};
    for (int f = 0; f < 8; ++f) { vector<double> v(nc); for (size_t i = 0; i < nc; ++i) v[i] = V.cands[i * 9 + f]; R->setv(names[f], std::move(v)); }
    const char* xn[3] = {"cand_a_voigt", "cand_s2_voigt", "cand_kind"};
    for (int f = 0; f < 3; ++f) { vector<double> v(nc); for (size_t i = 0; i < nc; ++i) v[i] = V.x[i * 3 + f]; R->setv(xn[f], std::move(v)); }
  }
  sbpublic::add_simple(data,o,*R);
  sbpublic::finish(*R,"voigt");
  R->finish_names();
  *out = reinterpret_cast<sb_result*>(guard.release());
  return SBC_OK;
}
}  

extern "C" {
SB_API void sb_opts_mass_init(sb_opts_mass* o) {
  if (!o) return;
  std::memset(o, 0, sizeof *o);
  o->struct_size = sizeof(sb_opts_mass); o->abi_version = SB_ABI_VERSION;
  o->threads = 0; o->check_all_windows = 1; o->seed = 1ULL; o->distribution_checks = 1;
  o->lp_iterations = SB_LP_ITER_FIRST; o->lp_retry_iterations = SB_LP_ITER_RETRY;
}
SB_API int sb_mass_gaussian(const sb_data* data, const sb_opts_mass* opts, sb_result** out) {
  return sbc::guarded([&] { return mass_gaussian_impl(data, opts, out); });
}
SB_API int sb_mass_voigt(const sb_data* data, const sb_opts_mass* opts, sb_result** out) {
  return sbc::guarded([&] { return mass_voigt_impl(data, opts, out, 2); });
}
SB_API int sb_test_famlev(const sb_data* data, const sb_opts_mass* opts, int family, sb_result** out) {
  if (family < 0 || family > 2) return sbc::fail(SBC_E_INPUT, "family must be 0, 1 or 2");
  return sbc::guarded([&] { return mass_voigt_impl(data, opts, out, family); });
}
SB_API int sb_test_ll_policy(const double* cands, size_t ncand, const char* internal, const double* summaries,
                             const sb_opts_mass* opts, sb_result** out) {
  return sbc::guarded([&] { return policy_impl(cands, ncand, internal, summaries, opts, out); });
}
}
