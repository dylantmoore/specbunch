// p2pilot.hpp: the Part 1 scale pilot, copied from cpp/part1/src/llcore.hpp (read-only source), 2026-11-04.
//   cf_recurrence      llcore.hpp lines 45-93   (the two omp-parallel pragmas removed: the fixed 64-chunk order is kept,
//                                                so sums are bitwise those of the parallel version at any thread count)
//   levrepmed_*, stable_band_filtered, ensemble_sigma   llcore.hpp lines 95-232 (unchanged)
//   ugrid_ctx, pilot_B llcore.hpp lines 234-269  (unchanged; used only for the group B-hat weights of the Y1 remedy)
// Everything sits in namespace p2pilot so it can be linked next to cpp/part2/friction.hpp.
// GCC 11 compatible (no <format>, no ranges).
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace p2pilot {
using std::vector;
static const double NaN = std::numeric_limits<double>::quiet_NaN();

static void cf_recurrence(const vector<double>& dZ, double t0, double Delta, int M, vector<double>& C, vector<double>& C2) {
  const size_t n = dZ.size();
  C.assign(M, 0.0); C2.assign(M, 0.0);
  constexpr int BLK = 512, NCH = 64;
  const size_t nblk = (n + BLK - 1) / BLK;
  const size_t per = (nblk + NCH - 1) / NCH;
  vector<double> P((size_t)NCH * M * 2, 0.0);
  {
    alignas(64) double cc[BLK], ss[BLK], cd[BLK], sd[BLK];
    for (int ch = 0; ch < NCH; ++ch) {
      double* sc = &P[(size_t)ch * M * 2]; double* sc2 = sc + M;
      const size_t b0 = ch * per, b1 = std::min(nblk, b0 + per);
      for (size_t b = b0; b < b1; ++b) {
        size_t i0 = b * BLK, i1 = std::min(n, i0 + BLK); int L = (int)(i1 - i0);
        for (int ii = 0; ii < L; ++ii) {
          double z = dZ[i0 + ii];
          cc[ii] = std::cos(t0 * z); ss[ii] = std::sin(t0 * z);
          cd[ii] = std::cos(Delta * z); sd[ii] = std::sin(Delta * z);
        }
        double a = 0, a2 = 0;
#pragma omp simd reduction(+:a,a2)
        for (int ii = 0; ii < L; ++ii) { a += cc[ii]; a2 += cc[ii] * cc[ii]; }
        sc[0] += a; sc2[0] += a2;
        for (int j = 1; j < M; ++j) {
          a = 0; a2 = 0;
#pragma omp simd reduction(+:a,a2)
          for (int ii = 0; ii < L; ++ii) {
            double c = cc[ii], s = ss[ii];
            double cn = c * cd[ii] - s * sd[ii];
            double sn = s * cd[ii] + c * sd[ii];
            cc[ii] = cn; ss[ii] = sn; a += cn; a2 += cn * cn;
          }
          sc[j] += a; sc2[j] += a2;
        }
      }
    }
  }
  for (int ch = 0; ch < NCH; ++ch)
    for (int j = 0; j < M; ++j) { C[j] += P[(size_t)ch * M * 2 + j]; C2[j] += P[(size_t)ch * M * 2 + M + j]; }
  const double inv_n = 1.0 / (double)n;
  for (int j = 0; j < M; ++j) { C[j] *= inv_n; C2[j] *= inv_n; }
}

static inline double range_step(double TL, double TU, int M) { return (TU - TL) / (M - 1); }

static double simple_median(vector<double> v) {
  if (v.empty()) return NaN;
  std::sort(v.begin(), v.end());
  size_t n = v.size();
  return (n % 2) ? v[(n - 1) / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

static double levrepmed_from_cached(const double* C, const double* t2, int M) {
  vector<double> gs(M);
  vector<double> betas; betas.reserve(M);
  for (int j = 1; j <= M - 5; ++j) {
    int n_partners = M - (j + 2);
    if (n_partners < 3) continue;
    double cj = C[j - 1], t2j = t2[j - 1];
    auto mg = [&](double beta) {
      int cnt = 0;
      for (int l = j + 3; l <= M; ++l) gs[cnt++] = C[l - 1] - cj * std::exp(-beta * (t2[l - 1] - t2j) * 0.5);
      std::sort(gs.begin(), gs.begin() + cnt);
      return (cnt % 2) ? gs[(cnt + 1) / 2 - 1] : 0.5 * (gs[cnt / 2 - 1] + gs[cnt / 2]);
    };
    double blo = 0.0001, bhi = 0.25;
    double mlo = mg(blo), mhi = mg(bhi);
    if (mlo * mhi > 0.0) continue;
    for (int it = 0; it < 30; ++it) {
      double bm = 0.5 * (blo + bhi);
      double mm = mg(bm);
      if (std::isnan(mm)) break;
      if (mlo * mm <= 0.0) { bhi = bm; mhi = mm; } else { blo = bm; mlo = mm; }
    }
    double bj = 0.5 * (blo + bhi);
    if (0.0 < bj && bj < 0.25) betas.push_back(bj);
  }
  if (betas.size() < 3) return NaN;
  return std::clamp(std::sqrt(simple_median(betas)), 0.005, 0.45);
}

static double levrepmed_on_band(const vector<double>& dZ, double TL, double TU, int M) {
  vector<double> C, C2;
  double D = range_step(TL, TU, M);
  cf_recurrence(dZ, TL, D, M, C, C2);
  vector<double> t2(M);
  for (int j = 0; j < M; ++j) { double t = (j == M - 1) ? TU : TL + j * D; t2[j] = t * t; }
  return levrepmed_from_cached(C.data(), t2.data(), M);
}

static inline double interp_cf(const vector<double>& Cd, double tlo, double dt_inv, int Md, double t) {
  double idx_f = (t - tlo) * dt_inv + 1.0;
  if (idx_f <= 1.0) return Cd[0];
  if (idx_f >= Md) return Cd[Md - 1];
  int lo = (int)std::floor(idx_f); lo = std::min(lo, Md - 1);
  double fr = idx_f - lo;
  return Cd[lo - 1] + fr * (Cd[lo] - Cd[lo - 1]);
}

static double stable_band_filtered(const vector<double>& dZ, double sig_wide) {
  const int n_tL = 40, n_tU = 5, M = 51, halfwin = 3;
  const double snr_c = 3.0, min_sf = 0.5;
  const double noise_floor = snr_c / std::sqrt((double)dZ.size());
  double sig0 = sig_wide;
  if (std::isnan(sig0)) return NaN;
  sig0 = std::clamp(sig0, 0.01, 0.45);
  double ts0 = std::sqrt(2.0) / sig0;
  double sig1 = levrepmed_on_band(dZ, 0.5 * ts0, 1.5 * ts0, 51);
  if (!std::isnan(sig1)) sig0 = std::clamp(sig1, 0.01, 0.45);
  double ts = std::sqrt(2.0) / sig0;
  double tdlo = std::max(0.3, 0.1 * ts), tdhi = 4.0 * ts;
  const int Md = 501;
  double dtd = (tdhi - tdlo) / (Md - 1);
  vector<double> Cd, C2d;
  cf_recurrence(dZ, tdlo, dtd, Md, Cd, C2d);
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

struct Pilot { double sig_filt, sig_iter0, sig_iter, sigma_p; };
static Pilot ensemble_sigma(const vector<double>& dZ) {
  Pilot P{};
  double wide = levrepmed_on_band(dZ, 2.0, 50.0, 81);
  P.sig_iter0 = wide;
  P.sig_filt = stable_band_filtered(dZ, wide);
  double si = wide;
  if (!std::isnan(si)) {
    si = std::clamp(si, 0.01, 0.45);
    for (int it = 0; it < 10; ++it) {
      double ts = std::sqrt(2.0) / si;
      double sn = levrepmed_on_band(dZ, 0.5 * ts, 1.5 * ts, 51);
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

constexpr double UGRID_LO = 0.20, UGRID_HI = 4.00, UGRID_DU = 0.05;
constexpr int UGRID_N = 77;
struct Ctx { size_t n; double sigma_p, t_start, Delta; vector<double> u, t, C, var, s; };
static Ctx ugrid_ctx(const vector<double>& dZ, double sp) {
  Ctx c; c.n = dZ.size(); c.sigma_p = sp;
  c.t_start = UGRID_LO / sp; c.Delta = UGRID_DU / sp;
  vector<double> C2;
  cf_recurrence(dZ, c.t_start, c.Delta, UGRID_N, c.C, C2);
  c.var.resize(UGRID_N); c.t.resize(UGRID_N); c.u.resize(UGRID_N); c.s.resize(UGRID_N);
  double tend = c.t_start + (UGRID_N - 1) * c.Delta;
  for (int j = 0; j < UGRID_N; ++j) {
    c.var[j] = C2[j] - c.C[j] * c.C[j];
    c.t[j] = (j == UGRID_N - 1) ? tend : c.t_start + j * ((tend - c.t_start) / (UGRID_N - 1));
    c.u[j] = (j == UGRID_N - 1) ? UGRID_HI : UGRID_LO + j * ((UGRID_HI - UGRID_LO) / (UGRID_N - 1));
    c.s[j] = std::sqrt(std::max(c.var[j] / (double)c.n, 1e-20));
  }
  return c;
}
static double pilot_B(const Ctx& c) {
  double num = 0, den = 0, inv_n = 1.0 / c.n;
  for (int j = 0; j < UGRID_N; ++j) {
    double uj = c.u[j];
    if (uj < 0.50 - 1e-12 || uj > 1.50 + 1e-12) continue;
    double sj2 = std::max(c.var[j] * inv_n, 1e-12), w = 1.0 / sj2, m = std::exp(-uj * uj / 2.0);
    num += w * m * c.C[j]; den += w * m * m;
  }
  if (!std::isfinite(num) || !std::isfinite(den) || den <= 0.0) return NaN;
  return std::max(1e-8, num / den);
}
}  // namespace p2pilot
