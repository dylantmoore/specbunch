// p1famv.hpp: the Voigt-friction analog of the headline Part 1 estimator, LLV.F3.S2, for the PUF run (2026-11-07). New code;
// llcore.hpp, p1w.hpp and every existing p1lib entry point are unchanged.
//
// What runs is cpp/part1/famlev (validated 2026-11-04 on the 600 robustness samples): the fits, the generic sup statistic and
// the runtime-width multiplier sums come from ../famfit.hpp, a byte copy of cpp/part1/famlev/src/famfit.hpp (assemble.py checks
// the bytes); famv_core below is famlev.cpp main()'s per-sample body for ONE family, transcribed line for line, with two
// additions only:
//   (1) probability weights: the pilot, the 77-point u-grid context, the extended-grid ECF and the multiplier sums use p1w's
//       weighted definitions (n -> n_r = round(n_eff), var_j -> n_r Var C_j, the multiplier perturbs C_j by
//       sum_i xi_i w_i (cos - C_j) / W). The extended-grid weighted multiplier sums are famfit's boot_sums_w with w_i
//       multiplying each row (boot_sums_ww below). w == nullptr runs famlev's unweighted code exactly.
//   (2) family 0 (G) skips the extended-grid scan, which G never uses (famlev computes it for C/V only).
// Families: 0 = G (LL.F3.S2; reproduces p1_spectral's LL.F3.S2 candidates bit for bit, used as a self-check), 1 = C, 2 = V
// (LLV.F3.S2: mu(t) = B exp(-a t - s2 t^2 / 2), a, s2 >= 0). Candidate tuple = llall/famlev layout
// [a, e, B_raw, sigma_eq, se_logB, se_log_sigma_eq, p, range_ok, null], plus x = [a_hat, s2_hat, kind] per candidate
// (kind: 0 interior, 1 boundary a = 0, 2 boundary s2 = 0, 3 interior via the KKT restart). Grid indices of the extended
// grid share the u-grid's (u = .20 + .05 j, up to u = 8).
#pragma once
#include "../famfit.hpp"

namespace p1v {

constexpr int MSCAN = 157;  // famlev: u in [0.2, 8]

// famfit boot_sums_w with probability weights: Y[b*MPd + j] = sum_i xi_{b,i} w_i cos(t_j dZ_i), S[b] = sum_i xi_{b,i} w_i.
// Same key, signs, groups of 64, chunking and accumulation order as famfit boot_sums_w / p1w boot_sums_w.
static void boot_sums_ww(const vector<double>& dZ, const vector<double>& w, double t0, double Delta, uint64_t seed, int M, int MPd,
                         vector<double>& Y, vector<double>& S) {
  const size_t n = dZ.size();
  const uint64_t key = splitmix64((uint64_t)seed + 7777777ULL);
  const size_t ng = (n + 63) / 64;
  constexpr int NCH = 48;
  const size_t per = (ng + NCH - 1) / NCH;
  const size_t YS = (size_t)BOOT_B * MPd;
  vector<double> Ych((size_t)NCH * YS, 0.0), Sch((size_t)NCH * BOOT_B, 0.0);
#pragma omp parallel
  {
    vector<double> Xv((size_t)64 * MPd);
    double wl[64];
#pragma omp for schedule(dynamic, 1)
    for (int ch = 0; ch < NCH; ++ch) {
      double* Yl = &Ych[(size_t)ch * YS]; double* Sl = &Sch[(size_t)ch * BOOT_B];
      const size_t g0 = ch * per, g1 = std::min(ng, g0 + per);
      for (size_t g = g0; g < g1; ++g) {
        size_t i0 = g * 64; int L = (int)std::min<size_t>(64, n - i0);
        for (int ii = 0; ii < 64; ++ii) {
          double* X = &Xv[(size_t)ii * MPd];
          if (ii >= L) { wl[ii] = 0.0; for (int j = 0; j < MPd; ++j) X[j] = 0.0; continue; }
          double z = dZ[i0 + ii], wi = w[i0 + ii]; wl[ii] = wi;
          double cl = std::cos(t0 * z), sl = std::sin(t0 * z), cd = std::cos(Delta * z), sd = std::sin(Delta * z);
          X[0] = wi * cl; double cp = cl, sp = sl;
          for (int j = 1; j < M; ++j) {
            double cn = cp * cd - sp * sd, sn = sp * cd + cp * sd;
            X[j] = wi * cn; cp = cn; sp = sn;
          }
          for (int j = M; j < MPd; ++j) X[j] = 0.0;
        }
        for (int b = 0; b < BOOT_B; ++b) {
          uint64_t wd = sign_word(key, b, g);
          double sacc = 0;
          for (int ii = 0; ii < L; ++ii) sacc += ((wd >> ii) & 1ULL) ? wl[ii] : -wl[ii];
          Sl[b] += sacc;
          double* yb = &Yl[(size_t)b * MPd];
          for (int jc = 0; jc < MPd; jc += 16) {
            double acc[16];
            for (int q = 0; q < 16; ++q) acc[q] = yb[jc + q];
            for (int ii = 0; ii < 64; ++ii) {
              double s = ((wd >> ii) & 1ULL) ? 1.0 : -1.0;
              const double* X = &Xv[(size_t)ii * MPd + jc];
#pragma omp simd
              for (int q = 0; q < 16; ++q) acc[q] += s * X[q];
            }
            for (int q = 0; q < 16; ++q) yb[jc + q] = acc[q];
          }
        }
      }
    }
  }
  Y.assign(YS, 0.0); S.assign(BOOT_B, 0.0);
  for (int ch = 0; ch < NCH; ++ch) {
    const double* Yl = &Ych[(size_t)ch * YS];
    for (size_t q = 0; q < YS; ++q) Y[q] += Yl[q];
    for (int b = 0; b < BOOT_B; ++b) S[b] += Sch[(size_t)ch * BOOT_B + b];
  }
}

struct Res {
  int internal = 0;            // 0 ok, 1 refused:pilot, 2 refused:no_valid_band
  double sigma_p = NaN, u_hi_ext = NaN;
  int NGE = 0, n_bands = 0, kkt_restarts = 0;
  double t_base = NaN, t_tab = NaN, t_boot = NaN, t_s2 = NaN, t_sel = NaN;
  vector<double> cands;        // n x 9 (famlev layout, NaN = null)
  vector<double> x;            // n x 3 [a_hat, s2_hat, kind] (G: [0, sigma^2, 1], as famlev)
};

// famlev.cpp main(), one sample, one family f. w == nullptr: famlev's unweighted path.
static void famv_core(const vector<double>& dZ, const vector<double>* w, uint64_t seed, int f, Res& O) {
  const double P0 = cpu_now();
  const size_t n0 = dZ.size();
  p1w::WS ws;
  if (w) ws = p1w::wsums(*w);
  Pilot P = w ? p1w::ensemble_sigma_w(dZ, *w, ws) : ensemble_sigma(dZ);
  O.sigma_p = P.sigma_p;
  if (!(std::isfinite(P.sigma_p) && P.sigma_p > 0)) { O.internal = 1; O.t_base = cpu_now() - P0; return; }
  // ---- G context (unchanged) and the extended C/V grid
  Ctx c = w ? p1w::ugrid_ctx_w(dZ, *w, ws, P.sigma_p) : ugrid_ctx(dZ, P.sigma_p);
  const size_t n = c.n;  // famlev: dZ.size(); weighted: n_r
  const double t0 = c.t[0], Dp = c.t[1] - c.t[0];
  int NGE = NG;
  vector<double> te, ve, sE, CE;
  if (f > 0) {
    vector<double> Ce, C2e, se(MSCAN), z(MSCAN);
    if (w) {
      p1w::cf_recurrence_w(dZ, *w, ws, t0, Dp, MSCAN, Ce, C2e, true);  // C2e = n_r Var C_j
      for (int j = 0; j < MSCAN; ++j) { se[j] = std::sqrt(std::max(C2e[j] / (double)n, 1e-20)); z[j] = Ce[j] / se[j]; }
    } else {
      cf_recurrence(dZ, t0, Dp, MSCAN, Ce, C2e);
      for (int j = 0; j < MSCAN; ++j) {
        se[j] = std::sqrt(std::max((C2e[j] - Ce[j] * Ce[j]) / (double)n, 1e-20));
        z[j] = Ce[j] / se[j];
      }
    }
    vector<double> zm = running_median5(z);
    int jc = -1;
    for (int j = 0; j < MSCAN; ++j) { if (zm[j] >= 2.0) jc = j; else break; }
    NGE = std::max(NG, jc + 1);
    te.resize(NGE); ve.resize(NGE); sE.resize(NGE); CE.assign(Ce.begin(), Ce.begin() + NGE);
    for (int j = 0; j < NGE; ++j) { te[j] = t0 + j * Dp; ve[j] = w ? C2e[j] : C2e[j] - Ce[j] * Ce[j]; sE[j] = se[j]; }
  }
  O.NGE = NGE; O.u_hi_ext = 0.20 + 0.05 * (NGE - 1);
  const int MPE = ((NGE + 15) / 16) * 16;
  vector<double> xs(NG);
  for (int j = 0; j < NG; ++j) xs[j] = -c.t[j] * c.t[j] / 2.0;
  O.t_base = cpu_now() - P0;
  // ---- bootstrap sums once (G: llcore boot_sums / p1w boot_sums_w; C/V on an extended grid: famfit width)
  double c0 = cpu_now();
  vector<double> YG, SG, YE, SE;
  if (f == 0 || NGE == NG) {
    if (w) p1w::boot_sums_w(dZ, *w, t0, Dp, seed, YG, SG); else boot_sums(dZ, t0, Dp, seed, YG, SG);
  } else {
    if (w) boot_sums_ww(dZ, *w, t0, Dp, seed, NGE, MPE, YE, SE); else boot_sums_w(dZ, t0, Dp, seed, NGE, MPE, YE, SE);
  }
  O.t_boot = cpu_now() - c0;
  const double inv_n = w ? 1.0 / ws.W : 1.0 / (double)n0;
  // family grid views
  const int NGf = (f == 0) ? NG : NGE;
  const double* Cf = (f == 0) ? c.C.data() : CE.data();
  const double* tf = (f == 0) ? c.t.data() : te.data();
  const double* vf = (f == 0) ? c.var.data() : ve.data();
  const double* sf = (f == 0) ? c.s.data() : sE.data();
  auto bid = [&](int a, int e) { return a * NGf + e; };
  // ---- table
  c0 = cpu_now();
  vector<BFX> T((size_t)NGf * NGf);
  vector<double> MU((size_t)NGf * NGf * NGf, 0.0);
  VDiag dg;
  for (int a = 0; a < NGf; ++a)
    for (int e = a + MINPTS - 1; e < NGf; ++e) {
      BFX o; int m = e - a + 1; bool v;
      if (f == 0) {
        v = fit_band<true>(2, c, xs.data(), c.C.data(), a, m, P.sigma_p, o.r);
        if (v) { o.a = 0.0; o.s2 = o.r.sigma * o.r.sigma; o.kind = 1; }
      } else {
        v = fit_fam<true>(f, Cf + a, tf + a, vf + a, m, n, P.sigma_p, o, &dg);
      }
      if (v) {
        if (f == 0) o.r.D = sup_std_resid(&c.C[a], &c.t[a], &c.s[a], m, o.r.B, o.r.sigma);
        else o.r.D = sup_std_resid_v(Cf + a, tf + a, sf + a, m, o.r.B, o.a, o.s2);
        double* mu = &MU[(size_t)bid(a, e) * NGf];
        for (int j = a; j <= e; ++j) mu[j] = (f == 0) ? mu_k0(o.r.B, o.r.sigma, c.t[j]) : mu_v(o.r.B, o.a, o.s2, tf[j]);
      }
      o.r.valid = v;
      T[bid(a, e)] = o;
    }
  O.kkt_restarts = dg.restarts;
  O.t_tab = cpu_now() - c0;
  // ---- bootstrap p
  c0 = cpu_now();
  vector<int> blist;
  for (int a = 0; a < NGf; ++a)
    for (int e = a + MINPTS - 1; e < NGf; ++e)
      if (T[bid(a, e)].r.valid && std::isfinite(T[bid(a, e)].r.D)) blist.push_back(bid(a, e));
  const int nbl = (int)blist.size();
  O.n_bands = nbl;
  vector<int> counts((size_t)NGf * NGf, 0);
#pragma omp parallel
  {
    vector<int> cl((size_t)NGf * NGf, 0);
    double Cs[GMAX];
#pragma omp for schedule(dynamic, 4)
    for (int b = 0; b < BOOT_B; ++b) {
      const bool wide = (f > 0 && NGE > NG);
      const double* yb = wide ? &YE[(size_t)b * MPE] : &YG[(size_t)b * MP];
      const double Sb = wide ? SE[b] : SG[b];
      for (int q = 0; q < nbl; ++q) {
        const int id = blist[q], a = id / NGf, e = id % NGf, m = e - a + 1;
        const double* mub = &MU[(size_t)id * NGf];
        for (int j = a; j <= e; ++j) Cs[j] = mub[j] + (yb[j] - Cf[j] * Sb) * inv_n;
        double Ds;
        if (f == 0) {
          BF r;
          if (!fit_band<false>(2, c, xs.data(), Cs, a, m, T[id].r.sigma, r)) continue;
          Ds = sup_std_resid(&Cs[a], &c.t[a], &c.s[a], m, r.B, r.sigma);
        } else {
          BFX o;
          if (!fit_fam<false>(f, Cs + a, tf + a, vf + a, m, n, T[id].r.sigma, o)) continue;
          Ds = sup_std_resid_v(Cs + a, tf + a, sf + a, m, o.r.B, o.a, o.s2);
        }
        if (Ds >= T[id].r.D) cl[id]++;
      }
    }
#pragma omp critical
    for (size_t q = 0; q < cl.size(); ++q) counts[q] += cl[q];
  }
  vector<double> PV((size_t)NGf * NGf, NaN);
  for (int id : blist) PV[id] = (1.0 + counts[id]) / (BOOT_B + 1.0);
  O.t_s2 = cpu_now() - c0;
  // ---- S2 ranking
  c0 = cpu_now();
  vector<int> hi, lo;
  for (int id : blist) (PV[id] >= 0.10 ? hi : lo).push_back(id);
  auto pts = [&](int id) { return id % NGf - id / NGf + 1; };
  std::stable_sort(hi.begin(), hi.end(), [&](int x, int y) {
    if (pts(x) != pts(y)) return pts(x) > pts(y);
    return x / NGf < y / NGf; });
  std::stable_sort(lo.begin(), lo.end(), [&](int x, int y) {
    if (PV[x] != PV[y]) return PV[x] > PV[y];
    if (pts(x) != pts(y)) return pts(x) > pts(y);
    return x / NGf < y / NGf; });
  vector<int> rk(hi); rk.insert(rk.end(), lo.begin(), lo.end());
  if (rk.empty()) O.internal = 2;
  O.cands.clear(); O.x.clear();
  O.cands.reserve(rk.size() * 9); O.x.reserve(rk.size() * 3);
  for (int id : rk) {
    const BF& r = T[id].r;
    const bool rok = std::isfinite(r.B) && r.B > 0.0 && r.B <= 1.5 && std::isfinite(r.sigma) && r.sigma >= 1e-3 && r.sigma <= 5.0;
    const double row[9] = {(double)(id / NGf), (double)(id % NGf), r.B, r.sigma, r.sevalid ? r.seB : NaN, r.sevalid ? r.seS : NaN,
                           PV[id], rok ? 1.0 : 0.0, NaN};
    O.cands.insert(O.cands.end(), row, row + 9);
    const double xr[3] = {T[id].a, T[id].s2, (double)T[id].kind};
    O.x.insert(O.x.end(), xr, xr + 3);
  }
  O.t_sel = cpu_now() - c0;
}

}  // namespace p1v
