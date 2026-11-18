// p1w.hpp: WEIGHTED variants of the llcore.hpp ECF machinery, for the PUF run (2026-11-06). New code; llcore.hpp is
// unchanged and still supplies the unweighted path, the fits and the selectors.
//
// Weighting (probability weights w_i > 0, W = sum w, W2 = sum w^2, Kish n_eff = W^2 / W2):
//   C_j    = sum_i w_i cos(t_j Z_i) / W                                   (weighted mean of cos; S likewise)
//   Var C_j = sum_i w_i^2 (cos(t_j Z_i) - C_j)^2 / W^2                     (exact variance of a weighted mean)
// Every llcore formula that uses (n, var_j) only through var_j / n keeps its meaning with n -> n_r = round(n_eff) and
// var_j -> n_r * Var C_j (Ctx.n, Ctx.var, Ctx.s, gqml_nr weights, fit_log weights, pilot_B, SNR1). The pilot noise
// floor 3 / sqrt(n) uses n_eff. The multiplier bootstrap perturbs C_j by sum_i xi_i w_i (cos - C_j) / W
// (unweighted: (1/n) sum_i xi_i (cos - C_j)), whose variance is Var C_j. All of it is scale-free in w.
// With w = 1 these reduce to the unweighted formulas (equal to rounding; p1lib keeps the exact llcore path for w = 1).
// Determinism: fixed chunks summed in chunk order, as llcore.hpp (bitwise identical at any thread count).
#pragma once
#include <utility>

namespace p1w {

struct WS { double W = 0, W2 = 0, neff = 0; size_t nr = 0; };
static WS wsums(const vector<double>& w) {
  WS s; long double a = 0, b = 0;
  for (double x : w) { a += x; b += (long double)x * x; }
  s.W = (double)a; s.W2 = (double)b; s.neff = (double)(a * a / b);
  s.nr = (size_t)std::max<long long>(1, std::llround(s.neff));
  return s;
}

// C[j] = weighted mean cos(t_j dZ); if want_var, Vn[j] = n_r * Var C_j. t_j = t0 + j Delta, j = 0..M-1.
static void cf_recurrence_w(const vector<double>& dZ, const vector<double>& w, const WS& ws, double t0, double Delta, int M,
                            vector<double>& C, vector<double>& Vn, bool want_var) {
  const size_t n = dZ.size();
  C.assign(M, 0.0); Vn.assign(M, 0.0);
  constexpr int BLK = 512, NCH = 64;
  const size_t nblk = (n + BLK - 1) / BLK;
  const size_t per = (nblk + NCH - 1) / NCH;
  const int NQ = want_var ? 3 : 1;
  vector<double> P((size_t)NCH * M * NQ, 0.0);
#pragma omp parallel
  {
    alignas(64) double cc[BLK], ss[BLK], cd[BLK], sd[BLK], ww[BLK], w2[BLK];
#pragma omp for schedule(dynamic, 1)
    for (int ch = 0; ch < NCH; ++ch) {
      double* s1 = &P[(size_t)ch * M * NQ]; double* q1 = s1 + M; double* q2 = s1 + 2 * M;
      const size_t b0 = ch * per, b1 = std::min(nblk, b0 + per);
      for (size_t b = b0; b < b1; ++b) {
        size_t i0 = b * BLK, i1 = std::min(n, i0 + BLK); int L = (int)(i1 - i0);
        for (int ii = 0; ii < L; ++ii) {
          double z = dZ[i0 + ii];
          cc[ii] = std::cos(t0 * z); ss[ii] = std::sin(t0 * z);
          cd[ii] = std::cos(Delta * z); sd[ii] = std::sin(Delta * z);
          ww[ii] = w[i0 + ii]; w2[ii] = ww[ii] * ww[ii];
        }
        for (int j = 0; j < M; ++j) {
          double a = 0, a1 = 0, a2 = 0;
          if (j > 0) {
            for (int ii = 0; ii < L; ++ii) {
              double c = cc[ii], s = ss[ii];
              cc[ii] = c * cd[ii] - s * sd[ii]; ss[ii] = s * cd[ii] + c * sd[ii];
            }
          }
          if (want_var) {
#pragma omp simd reduction(+:a,a1,a2)
            for (int ii = 0; ii < L; ++ii) { double c = cc[ii]; a += ww[ii] * c; a1 += w2[ii] * c; a2 += w2[ii] * c * c; }
            s1[j] += a; q1[j] += a1; q2[j] += a2;
          } else {
#pragma omp simd reduction(+:a)
            for (int ii = 0; ii < L; ++ii) a += ww[ii] * cc[ii];
            s1[j] += a;
          }
        }
      }
    }
  }
  vector<double> S1(M, 0.0), Q1(M, 0.0), Q2(M, 0.0);
  for (int ch = 0; ch < NCH; ++ch) {
    const double* s1 = &P[(size_t)ch * M * NQ];
    for (int j = 0; j < M; ++j) {
      S1[j] += s1[j];
      if (want_var) { Q1[j] += s1[M + j]; Q2[j] += s1[2 * M + j]; }
    }
  }
  const double iW = 1.0 / ws.W;
  for (int j = 0; j < M; ++j) {
    C[j] = S1[j] * iW;
    if (want_var) {
      double v = (Q2[j] - 2.0 * C[j] * Q1[j] + C[j] * C[j] * ws.W2) * iW * iW;
      Vn[j] = std::max(v, 0.0) * (double)ws.nr;
    }
  }
}

// ---------------------------------------------------------------- pilot (llcore lines 136-232, weighted ECF)
static double levrepmed_on_band_w(const vector<double>& dZ, const vector<double>& w, const WS& ws, double TL, double TU, int M) {
  vector<double> C, Vn;
  double D = range_step(TL, TU, M);
  cf_recurrence_w(dZ, w, ws, TL, D, M, C, Vn, false);
  vector<double> t2(M);
  for (int j = 0; j < M; ++j) { double t = (j == M - 1) ? TU : TL + j * D; t2[j] = t * t; }
  return levrepmed_from_cached(C.data(), t2.data(), M);
}

// llcore stable_band_filtered with the weighted ECF; noise floor 3 / sqrt(n_eff)
static double stable_band_filtered_w(const vector<double>& dZ, const vector<double>& w, const WS& ws, double sig_wide) {
  const int n_tL = 40, n_tU = 5, M = 51, halfwin = 3;
  const double snr_c = 3.0, min_sf = 0.5;
  const double noise_floor = snr_c / std::sqrt(ws.neff);
  double sig0 = sig_wide;
  if (std::isnan(sig0)) return NaN;
  sig0 = std::clamp(sig0, 0.01, 0.45);
  double ts0 = std::sqrt(2.0) / sig0;
  double sig1 = levrepmed_on_band_w(dZ, w, ws, 0.5 * ts0, 1.5 * ts0, 51);
  if (!std::isnan(sig1)) sig0 = std::clamp(sig1, 0.01, 0.45);
  double ts = std::sqrt(2.0) / sig0;
  double tdlo = std::max(0.3, 0.1 * ts), tdhi = 4.0 * ts;
  const int Md = 501;
  double dtd = (tdhi - tdlo) / (Md - 1);
  vector<double> Cd, Vd;
  cf_recurrence_w(dZ, w, ws, tdlo, dtd, Md, Cd, Vd, false);
  double dt_inv = 1.0 / dtd;
  double tLlo = std::max(0.5, 0.15 * ts), tLhi = 1.1 * ts, tUlo = 1.5 * ts, tUhi = 3.0 * ts;
  auto rng = [](double a, double b, int L, int i) { return (i == L - 1) ? b : a + i * ((b - a) / (L - 1)); };
  double best_sig = NaN, best_score = INFINITY;
  vector<double> sigv(n_tL), sfrac(n_tL), Cb(M), t2b(M);
  for (int iu = 0; iu < n_tU; ++iu) {
    double tU = rng(tUlo, tUhi, n_tU, iu);
    for (int il = 0; il < n_tL; ++il) {
      double tL = rng(tLlo, tLhi, n_tL, il);
      if (tL >= tU - 1.0) { sigv[il] = NaN; sfrac[il] = 0.0; continue; }
      double dtb = (tU - tL) / (M - 1);
      int nab = 0;
      for (int m = 0; m < M; ++m) {
        double tm = tL + m * dtb;
        Cb[m] = interp_cf(Cd, tdlo, dt_inv, Md, tm);
        t2b[m] = tm * tm;
        if (std::fabs(Cb[m]) > noise_floor) ++nab;
      }
      sfrac[il] = (double)nab / M;
      sigv[il] = levrepmed_from_cached(Cb.data(), t2b.data(), M);
    }
    for (int il = halfwin; il < n_tL - halfwin; ++il) {
      if (sfrac[il] < min_sf) continue;
      bool ok = true; double lo = INFINITY, hi = -INFINITY;
      for (int jj = il - halfwin; jj <= il + halfwin; ++jj) {
        double sv = sigv[jj];
        if (std::isnan(sv)) { ok = false; break; }
        lo = std::min(lo, sv); hi = std::max(hi, sv);
      }
      if (!ok) continue;
      double score = (hi - lo) / std::max(sfrac[il], 0.01);
      if (score < best_score) { best_score = score; best_sig = sigv[il]; }
    }
  }
  return std::isnan(best_sig) ? sig0 : std::clamp(best_sig, 0.005, 0.45);
}

static Pilot ensemble_sigma_w(const vector<double>& dZ, const vector<double>& w, const WS& ws) {
  Pilot P{};
  double wide = levrepmed_on_band_w(dZ, w, ws, 2.0, 50.0, 81);
  P.sig_iter0 = wide;
  P.sig_filt = stable_band_filtered_w(dZ, w, ws, wide);
  double si = wide;
  if (!std::isnan(si)) {
    si = std::clamp(si, 0.01, 0.45);
    for (int it = 0; it < 10; ++it) {
      double ts = std::sqrt(2.0) / si;
      double sn = levrepmed_on_band_w(dZ, w, ws, 0.5 * ts, 1.5 * ts, 51);
      if (std::isnan(sn)) break;
      if (std::fabs(sn - si) / std::max(si, 0.01) < 0.01) { si = sn; break; }
      si = sn;
    }
  }
  P.sig_iter = si;
  if (std::isnan(P.sig_filt) && std::isnan(si)) P.sigma_p = NaN;
  else if (std::isnan(P.sig_filt)) P.sigma_p = si;
  else if (std::isnan(si)) P.sigma_p = P.sig_filt;
  else P.sigma_p = 0.5 * (P.sig_filt + si);
  return P;
}

// llcore ugrid_ctx with the weighted ECF: c.n = n_r, c.var = n_r Var C_j
static Ctx ugrid_ctx_w(const vector<double>& dZ, const vector<double>& w, const WS& ws, double sp) {
  Ctx c; c.n = ws.nr; c.sigma_p = sp;
  c.t_start = UGRID_LO / sp; c.Delta = UGRID_DU / sp;
  cf_recurrence_w(dZ, w, ws, c.t_start, c.Delta, UGRID_N, c.C, c.var, true);
  c.t.resize(UGRID_N); c.u.resize(UGRID_N); c.s.resize(UGRID_N);
  double tend = c.t_start + (UGRID_N - 1) * c.Delta;
  for (int j = 0; j < UGRID_N; ++j) {
    c.t[j] = (j == UGRID_N - 1) ? tend : c.t_start + j * ((tend - c.t_start) / (UGRID_N - 1));
    c.u[j] = (j == UGRID_N - 1) ? UGRID_HI : UGRID_LO + j * ((UGRID_HI - UGRID_LO) / (UGRID_N - 1));
    c.s[j] = std::sqrt(std::max(c.var[j] / (double)c.n, 1e-20));
  }
  return c;
}

// llcore boot_sums with weights: Y[b*MP + j] = sum_i xi_{b,i} w_i cos(t_j dZ_i), S[b] = sum_i xi_{b,i} w_i.
// Same Rademacher stream (sign_word) as the unweighted path; fixed chunks summed in chunk order.
static void boot_sums_w(const vector<double>& dZ, const vector<double>& w, double t0, double Delta, uint64_t seed,
                        vector<double>& Y, vector<double>& S) {
  const size_t n = dZ.size();
  const uint64_t key = splitmix64((uint64_t)seed + 7777777ULL);
  const size_t ng = (n + 63) / 64;
  constexpr int NCH = 48;
  const size_t per = (ng + NCH - 1) / NCH;
  const size_t YS = (size_t)BOOT_B * MP;
  vector<double> Ych((size_t)NCH * YS, 0.0), Sch((size_t)NCH * BOOT_B, 0.0);
#pragma omp parallel
  {
    alignas(64) double X[64][MP];
    double wl[64];
#pragma omp for schedule(dynamic, 1)
    for (int ch = 0; ch < NCH; ++ch) {
      double* Yl = &Ych[(size_t)ch * YS]; double* Sl = &Sch[(size_t)ch * BOOT_B];
      const size_t g0 = ch * per, g1 = std::min(ng, g0 + per);
      for (size_t g = g0; g < g1; ++g) {
        size_t i0 = g * 64; int L = (int)std::min<size_t>(64, n - i0);
        for (int ii = 0; ii < 64; ++ii) {
          if (ii >= L) { wl[ii] = 0.0; for (int j = 0; j < MP; ++j) X[ii][j] = 0.0; continue; }
          double z = dZ[i0 + ii], wi = w[i0 + ii]; wl[ii] = wi;
          double cl = std::cos(t0 * z), sl = std::sin(t0 * z), cd = std::cos(Delta * z), sd = std::sin(Delta * z);
          X[ii][0] = wi * cl; double cp = cl, sp = sl;
          for (int j = 1; j < UGRID_N; ++j) {
            double cn = cp * cd - sp * sd, sn = sp * cd + cp * sd;
            X[ii][j] = wi * cn; cp = cn; sp = sn;
          }
          for (int j = UGRID_N; j < MP; ++j) X[ii][j] = 0.0;
        }
        for (int b = 0; b < BOOT_B; ++b) {
          uint64_t wd = sign_word(key, b, g);
          double sacc = 0;
          for (int ii = 0; ii < L; ++ii) sacc += ((wd >> ii) & 1ULL) ? wl[ii] : -wl[ii];
          Sl[b] += sacc;
          double* yb = &Yl[(size_t)b * MP];
          for (int jc = 0; jc < MP; jc += 16) {
            double acc[16];
            for (int q = 0; q < 16; ++q) acc[q] = yb[jc + q];
            for (int ii = 0; ii < 64; ++ii) {
              double s = ((wd >> ii) & 1ULL) ? 1.0 : -1.0;
#pragma omp simd
              for (int q = 0; q < 16; ++q) acc[q] += s * X[ii][jc + q];
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

}  // namespace p1w
