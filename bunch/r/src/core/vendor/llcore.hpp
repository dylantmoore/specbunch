// llcore.hpp: deterministic copy of the verified gqcsn.cpp machinery (cppest/spec/gqcsn.cpp lines 13-435,
// 490-497) and the llwin.cpp band fits/selectors (cppest/llwin/llwin.cpp lines 25-186), for RESTRICT_ALL.
// Changes vs the originals (results otherwise identical up to floating-point summation order):
//   cf_recurrence, boot_sums: per-thread partials + omp critical (order-dependent) replaced by a fixed
//   number of chunks summed in chunk order, so results are bitwise identical at any thread count.
#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <omp.h>

using std::vector;
static const double NaN = std::numeric_limits<double>::quiet_NaN();
static double now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// ---------------------------------------------------------------- npy
static vector<double> read_npy(const std::string& path) {
  std::ifstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("cannot open " + path);
  char magic[6]; f.read(magic, 6);
  if (std::memcmp(magic + 1, "NUMPY", 5) != 0) throw std::runtime_error("not npy");
  unsigned char ver[2]; f.read((char*)ver, 2);
  uint32_t hl = 0;
  if (ver[0] == 1) { uint16_t h; f.read((char*)&h, 2); hl = h; } else { f.read((char*)&hl, 4); }
  std::string hdr(hl, ' '); f.read(hdr.data(), hl);
  if (hdr.find("<f8") == std::string::npos || hdr.find("False") == std::string::npos)
    throw std::runtime_error("npy must be <f8 C-order");
  auto p0 = f.tellg(); f.seekg(0, std::ios::end); auto p1 = f.tellg(); f.seekg(p0);
  size_t n = (size_t)(p1 - p0) / 8;
  vector<double> z(n); f.read((char*)z.data(), n * 8);
  return z;
}

// ---------------------------------------------------------------- CF via recurrence
// C[j] = mean cos(t_j dZ), C2[j] = mean cos^2, t_j = t0 + j*Delta, j=0..M-1.
static void cf_recurrence(const vector<double>& dZ, double t0, double Delta, int M,
                          vector<double>& C, vector<double>& C2) {
  // Deterministic: fixed NCH chunks of 512-obs blocks (independent of thread count); each chunk
  // sums its blocks in order into its own slot; slots are added in chunk order.
  const size_t n = dZ.size();
  C.assign(M, 0.0); C2.assign(M, 0.0);
  constexpr int BLK = 512, NCH = 64;
  const size_t nblk = (n + BLK - 1) / BLK;
  const size_t per = (nblk + NCH - 1) / NCH;
  vector<double> P((size_t)NCH * M * 2, 0.0);
#pragma omp parallel
  {
    alignas(64) double cc[BLK], ss[BLK], cd[BLK], sd[BLK];
#pragma omp for schedule(dynamic, 1)
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

// Julia range(TL, TU, length=M) first step
static inline double range_step(double TL, double TU, int M) { return (TU - TL) / (M - 1); }

// ---------------------------------------------------------------- LevRepMed
static double simple_median(vector<double> v) {
  if (v.empty()) return NaN;
  std::sort(v.begin(), v.end());
  size_t n = v.size();
  return (n % 2) ? v[(n - 1) / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// _levrepmed_from_cached (identical algorithm to _levrepmed_on_band given C)
static double levrepmed_from_cached(const double* C, const double* t2, int M) {
  vector<double> gs(M);
  vector<double> betas; betas.reserve(M);
  for (int j = 1; j <= M - 5; ++j) {  // Julia 1-based j
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

// _stable_band_filtered with defaults n_tL=40,n_tU=5,M=51,halfwin=3,snr_c=3,min_signal_frac=.5
// sig_wide = _levrepmed_on_band(Z,k,2,50,81) (precomputed, deterministic)
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

// ---------------------------------------------------------------- u-grid context
constexpr double UGRID_LO = 0.20, UGRID_HI = 4.00, UGRID_DU = 0.05;
constexpr int UGRID_N = 77, MP = 80;  // MP = padded row length
constexpr int BOOT_B = 499;

struct Ctx {
  size_t n; double sigma_p, t_start, Delta;
  vector<double> u, t, C, var, s;
};
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

static vector<double> running_median5(const vector<double>& x) {
  int n = (int)x.size(); vector<double> out(n);
  for (int j = 0; j < n; ++j) {
    int lo = std::max(0, j - 2), hi = std::min(n - 1, j + 2);
    vector<double> b(x.begin() + lo, x.begin() + hi + 1);
    out[j] = simple_median(b);
  }
  return out;
}

enum Policy { Fixed2p5 = 0, SNR1 = 1, Info99 = 2 };
static const char* POLNAME[3] = {"Fixed2p5", "SNR1", "Info99"};

// returns Julia 1-based index, or -1
static int choose_uH_idx(Policy p, const Ctx& c, double Bp, const vector<double>& ssm) {
  const int M = UGRID_N;
  if (p == Fixed2p5) {
    int idx = -1;
    for (int j = 0; j < M; ++j) if (c.u[j] <= 2.50 + 1e-12) idx = j + 1;
    return idx;
  } else if (p == SNR1) {
    int best = -1;
    for (int j = 0; j < M; ++j) {
      double mu = Bp * std::exp(-c.u[j] * c.u[j] / 2.0), q = ssm[j];
      if (q <= 0.0 || !std::isfinite(q)) continue;
      double sa = mu / q, sp = c.u[j] * c.u[j] * mu / q;
      if (std::max(sa, sp) >= 1.0) best = j + 1;
    }
    return best;
  } else {
    double inv_n = 1.0 / c.n, total = 0; vector<double> I(M);
    for (int j = 0; j < M; ++j) {
      double sj2 = std::max(c.var[j] * inv_n, 1e-12), w = 1.0 / sj2;
      double mu = Bp * std::exp(-c.u[j] * c.u[j] / 2.0);
      double a = c.u[j] * c.u[j] * mu;
      I[j] = w * (mu * mu + a * a); total += I[j];
    }
    if (!std::isfinite(total) || total <= 0.0) return -1;
    double cum = 0;
    for (int j = 0; j < M; ++j) { cum += I[j]; if (cum / total >= 0.99) return j + 1; }
    return M;
  }
}

static vector<int> candidate_lo(int uH, const vector<double>& u) {
  double cap = u[uH - 1] - 0.50; vector<int> c;
  for (int j = 5; j <= uH; j += 2) {
    if (u[j - 1] > cap + 1e-12) break;
    if (uH - j + 1 >= 11) c.push_back(j);
  }
  return c;
}

// ---------------------------------------------------------------- GQML NR (K=0)
struct NR { double B, sigma; bool conv; };
static NR gqml_nr(const double* C, const double* t, const double* var, int M, size_t n,
                  double sigma_init, int max_iter = 100, double tol = 1e-8) {
  double w[UGRID_N];
  for (int j = 0; j < M; ++j) w[j] = var[j] > 1e-12 ? (double)n / var[j] : (double)n;
  double sxx = 0, sx = 0, sy = 0, sxy = 0;
  for (int j = 0; j < M; ++j) {
    double x = -t[j] * t[j] / 2.0, y = std::log(std::max(C[j], 1e-10));
    sx += x; sy += y; sxx += x * x; sxy += x * y;
  }
  double Md = (double)M, denom = Md * sxx - sx * sx, alpha, s2i;
  if (std::fabs(denom) > 1e-30) {
    double b2 = (Md * sxy - sx * sy) / denom;
    alpha = (sy - b2 * sx) / Md;
    s2i = std::max((Md * sxy - sx * sy) / denom, sigma_init * sigma_init);
  } else { alpha = std::log(std::max(C[0], 1e-10)); s2i = sigma_init * sigma_init; }
  double psi = 0.5 * std::log(s2i);
  bool conv = false;
  for (int it = 0; it < max_iter; ++it) {
    double s2 = std::exp(2.0 * psi);
    double ga = 0, gp = 0, Haa = 0, Hap = 0, Hpp = 0, Ga = 0, Gap = 0, Gpp = 0;
    for (int j = 0; j < M; ++j) {
      double t2 = t[j] * t[j], aj = s2 * t2;
      double eta = std::clamp(alpha - aj / 2.0, -500.0, 500.0), mj = std::exp(eta);
      double rj = C[j] - mj, wj = w[j], wm = wj * mj, wm2 = wm * mj;
      ga += wm * rj; gp -= wm * aj * rj;
      Haa += wm * (rj - mj); Hap += wm * aj * (mj - rj);
      Hpp += wj * aj * mj * ((aj - 2.0) * rj - aj * mj);
      Ga -= wm2; Gap += wm2 * aj; Gpp -= wm2 * aj * aj;
    }
    double detf = Haa * Hpp - Hap * Hap, h11, h12, h22, det;
    if (Haa < 0.0 && detf > 1e-30) { h11 = Haa; h12 = Hap; h22 = Hpp; det = detf; }
    else { h11 = Ga; h12 = Gap; h22 = Gpp; det = Ga * Gpp - Gap * Gap; }
    if (std::fabs(det) < 1e-30) break;
    double id = 1.0 / det;
    double da = id * (h22 * ga - h12 * gp), dp = id * (-h12 * ga + h11 * gp);
    if (!std::isfinite(da) || !std::isfinite(dp)) break;
    alpha -= da; psi -= dp;
    if (da * da + dp * dp < tol * tol * (1.0 + alpha * alpha + psi * psi)) { conv = true; break; }
  }
  return {std::exp(std::clamp(alpha, -500.0, 500.0)), std::exp(psi), conv};
}

static inline double mu_k0(double B, double s, double t) { return B * std::exp(-s * s * t * t * 0.5); }
static double sup_std_resid(const double* C, const double* t, const double* s, int M, double B, double sg) {
  double D = 0;
  for (int j = 0; j < M; ++j) { double d = std::fabs(C[j] - mu_k0(B, sg, t[j])) / s[j]; if (d > D) D = d; }
  return D;
}

// ---------------------------------------------------------------- bootstrap multiplier sums
static inline uint64_t splitmix64(uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}
// 64 Rademacher signs for (draw b, obs group g): bit=1 -> +1, bit=0 -> -1
static inline uint64_t sign_word(uint64_t key, int b, uint64_t g) {
  return splitmix64(key ^ splitmix64(((uint64_t)b << 40) ^ g));
}

// Y[b*MP + j] = sum_i xi_{b,i} cos(t_j dZ_i), S[b] = sum_i xi_{b,i}
// Deterministic: groups of 64 obs are split into NCH fixed chunks (independent of thread count);
// each chunk accumulates in group order into its own slot; slots are added in chunk order.
static void boot_sums(const vector<double>& dZ, double t0, double Delta, uint64_t seed,
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
#pragma omp for schedule(dynamic, 1)
    for (int ch = 0; ch < NCH; ++ch) {
      double* Yl = &Ych[(size_t)ch * YS]; double* Sl = &Sch[(size_t)ch * BOOT_B];
      const size_t g0 = ch * per, g1 = std::min(ng, g0 + per);
      for (size_t g = g0; g < g1; ++g) {
        size_t i0 = g * 64; int L = (int)std::min<size_t>(64, n - i0);
        for (int ii = 0; ii < 64; ++ii) {
          if (ii >= L) { for (int j = 0; j < MP; ++j) X[ii][j] = 0.0; continue; }
          double z = dZ[i0 + ii];
          double cl = std::cos(t0 * z), sl = std::sin(t0 * z), cd = std::cos(Delta * z), sd = std::sin(Delta * z);
          X[ii][0] = cl; double cp = cl, sp = sl;
          for (int j = 1; j < UGRID_N; ++j) {
            double cn = cp * cd - sp * sd, sn = sp * cd + cp * sd;
            X[ii][j] = cn; cp = cn; sp = sn;
          }
          for (int j = UGRID_N; j < MP; ++j) X[ii][j] = 0.0;
        }
        const uint64_t mask = (L == 64) ? ~0ULL : ((1ULL << L) - 1);
        for (int b = 0; b < BOOT_B; ++b) {
          uint64_t w = sign_word(key, b, g);
          Sl[b] += 2.0 * __builtin_popcountll(w & mask) - (double)L;
          double* yb = &Yl[(size_t)b * MP];
          for (int jc = 0; jc < MP; jc += 16) {
            double acc[16];
            for (int q = 0; q < 16; ++q) acc[q] = yb[jc + q];
            for (int ii = 0; ii < 64; ++ii) {
              double s = ((w >> ii) & 1ULL) ? 1.0 : -1.0;
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

static std::string jnum(double x) {
  if (!std::isfinite(x)) return "null";
  char b[40]; std::snprintf(b, sizeof b, "%.17g", x); return b;
}
static std::string fnum(double x) {
  if (!std::isfinite(x)) return "";
  char b[40]; std::snprintf(b, sizeof b, "%.17g", x); return b;
}
static double cpu_now() {
  rusage r; getrusage(RUSAGE_SELF, &r);
  return r.ru_utime.tv_sec + 1e-6 * r.ru_utime.tv_usec + r.ru_stime.tv_sec + 1e-6 * r.ru_stime.tv_usec;
}

constexpr int NG = UGRID_N;   // 77
constexpr int MINPTS = 11;
static inline int bidx(int a, int e) { return a * NG + e; }

// ---------------------------------------------------------------- band fits
struct BF {
  double B = NaN, sigma = NaN, logB = NaN, logS = NaN, seB = NaN, seS = NaN, D = NaN;
  bool valid = false;   // fit-level validity (spec)
  bool sevalid = false; // finite positive standard errors
};

// F1/F2: (weighted) least squares of log C on (1, x), x = -t^2/2. Returns validity.
// SE: F1 sandwich with v_j = Var_j/(n C_j^2); F2 model-based (X'WX)^{-1}, w_j = 1/v_j.
template <bool WEIGHTED, bool WANT_SE>
static inline bool fit_log(const double* C, const double* x, const double* var, int m, double nd, BF& r) {
  double y[NG], w[NG];
  for (int j = 0; j < m; ++j) {
    if (!(C[j] > 0.0)) return false;
    y[j] = std::log(C[j]);
  }
  if (WEIGHTED || WANT_SE)
    for (int j = 0; j < m; ++j) w[j] = nd * C[j] * C[j] / std::max(var[j], 1e-12);  // 1/v_j
  double W = 0, sx = 0, sy = 0;
  if (WEIGHTED) {
    for (int j = 0; j < m; ++j) { W += w[j]; sx += w[j] * x[j]; sy += w[j] * y[j]; }
  } else {
    W = m;
    for (int j = 0; j < m; ++j) { sx += x[j]; sy += y[j]; }
  }
  const double xb = sx / W, yb = sy / W;
  double Sxx = 0, Sxy = 0;
  if (WEIGHTED) {
    for (int j = 0; j < m; ++j) { double dx = x[j] - xb; Sxx += w[j] * dx * dx; Sxy += w[j] * dx * (y[j] - yb); }
  } else {
    for (int j = 0; j < m; ++j) { double dx = x[j] - xb; Sxx += dx * dx; Sxy += dx * (y[j] - yb); }
  }
  if (!(Sxx > 0.0)) return false;
  const double b = Sxy / Sxx, a = yb - b * xb;
  if (!(b > 0.0) || !std::isfinite(a) || !std::isfinite(b)) return false;
  r.logB = a; r.B = std::exp(a); r.sigma = std::sqrt(b); r.logS = 0.5 * std::log(b);
  if (!std::isfinite(r.B) || !std::isfinite(r.sigma)) return false;
  if (WANT_SE) {
    double vB, vb;
    if (WEIGHTED) { vb = 1.0 / Sxx; vB = 1.0 / W + xb * xb / Sxx; }
    else {
      // b_hat = sum c_j y_j, c_j = (x_j - xb)/Sxx; a_hat = sum d_j y_j, d_j = 1/m - xb c_j
      vb = 0; vB = 0;
      for (int j = 0; j < m; ++j) {
        double cj = (x[j] - xb) / Sxx, dj = 1.0 / m - xb * cj, vj = 1.0 / w[j];
        vb += cj * cj * vj; vB += dj * dj * vj;
      }
    }
    r.seB = std::sqrt(vB); r.seS = std::sqrt(vb) / (2.0 * b);
    r.sevalid = std::isfinite(r.seB) && std::isfinite(r.seS) && r.seB > 0 && r.seS > 0;
  }
  return true;
}

// F3: gqml_nr (unchanged), SE from Gauss-Newton weighted information in (log B, log sigma).
template <bool WANT_SE>
static inline bool fit_lev(const double* C, const double* t, const double* var, int m, size_t n,
                           double sigma_init, BF& r) {
  NR q = gqml_nr(C, t, var, m, n, sigma_init);
  if (!(q.conv && std::isfinite(q.B) && std::isfinite(q.sigma))) return false;
  r.B = q.B; r.sigma = q.sigma; r.logB = std::log(q.B); r.logS = std::log(q.sigma);
  if (WANT_SE) {
    const double s2 = q.sigma * q.sigma;
    double J11 = 0, J12 = 0, J22 = 0;
    for (int j = 0; j < m; ++j) {
      double wj = var[j] > 1e-12 ? (double)n / var[j] : (double)n;
      double aj = s2 * t[j] * t[j], mj = q.B * std::exp(-aj * 0.5), wm2 = wj * mj * mj;
      J11 += wm2; J12 -= wm2 * aj; J22 += wm2 * aj * aj;
    }
    double det = J11 * J22 - J12 * J12;
    r.seB = std::sqrt(J22 / det); r.seS = std::sqrt(J11 / det);
    r.sevalid = det > 0 && std::isfinite(r.seB) && std::isfinite(r.seS) && r.seB > 0 && r.seS > 0;
  }
  return true;
}

// Generic dispatch. F: 0=F1, 1=F2, 2=F3.
template <bool WANT_SE>
static inline bool fit_band(int F, const Ctx& c, const double* xs, const double* Cv, int a, int m,
                            double sigma_init, BF& r) {
  if (F == 0) return fit_log<false, WANT_SE>(Cv + a, xs + a, c.var.data() + a, m, (double)c.n, r);
  if (F == 1) return fit_log<true, WANT_SE>(Cv + a, xs + a, c.var.data() + a, m, (double)c.n, r);
  return fit_lev<WANT_SE>(Cv + a, c.t.data() + a, c.var.data() + a, m, c.n, sigma_init, r);
}

static const char* FNAME[3] = {"F1", "F2", "F3"};
static const char* FDESC[3] = {"log OLS", "log WLS", "levels WLS (gqml_nr)"};

// ---------------------------------------------------------------- selection result
struct Pick {
  std::string status = "refused:no_valid_band";
  double B = NaN, sigma = NaN;
  int a = -1, e = -1;
  double extra = NaN;  // p (S2) or plateau score (S3)
};
static void finalize(Pick& p, const BF& r) {
  bool ok = std::isfinite(r.B) && r.B > 0.0 && r.B <= 1.5 && std::isfinite(r.sigma) && r.sigma >= 1e-3 && r.sigma <= 5.0;
  if (!ok) { p.status = "refused:invalid_fit"; return; }
  p.status = "ok"; p.B = std::clamp(r.B, 0.0, 1.0); p.sigma = r.sigma;
}

// S1 Lepski over u_L in [0, e-10] (0-based), fixed u_H index e. Global rule: lowest valid u_L
// consistent with every narrower valid band. Also returns the sequential (stop-at-first-failure)
// choice for diagnostics (aseq).
static Pick lepski(const vector<BF>& T, int e, double kappa, int& aseq) {
  Pick p; aseq = -1;
  const int amax = e - (MINPTS - 1);
  auto usable = [&](int a) { const BF& r = T[bidx(a, e)]; return r.valid && r.sevalid; };
  auto consistent = [&](int a) {
    const BF& r = T[bidx(a, e)];
    for (int a2 = a + 1; a2 <= amax; ++a2) {
      if (!usable(a2)) continue;
      const BF& q = T[bidx(a2, e)];
      if (!(std::fabs(r.logB - q.logB) <= kappa * q.seB)) return false;
      if (!(std::fabs(r.logS - q.logS) <= kappa * q.seS)) return false;
    }
    return true;
  };
  int best = -1;
  for (int a = amax; a >= 0; --a) if (usable(a) && consistent(a)) best = a;  // keep lowest
  // sequential diagnostic
  for (int a = amax; a >= 0; --a) {
    if (!usable(a)) continue;
    if (consistent(a)) aseq = a; else break;
  }
  if (best < 0) return p;
  p.a = best; p.e = e;
  finalize(p, T[bidx(best, e)]);
  return p;
}

// S3 plateau: windows of 7 grid-consecutive u_L candidates, all valid.
static Pick plateau(const vector<BF>& T, int e) {
  Pick p; p.status = "refused:no_window";
  const int amax = e - (MINPTS - 1);
  double best = INFINITY; int bi = -1;
  for (int i = 0; i + 6 <= amax; ++i) {
    bool ok = true; double lb = INFINITY, hb = -INFINITY, ls = INFINITY, hs = -INFINITY;
    for (int a = i; a <= i + 6; ++a) {
      const BF& r = T[bidx(a, e)];
      if (!r.valid) { ok = false; break; }
      lb = std::min(lb, r.logB); hb = std::max(hb, r.logB); ls = std::min(ls, r.logS); hs = std::max(hs, r.logS);
    }
    if (!ok) continue;
    double sc = (hb - lb) + (hs - ls);
    if (sc < best) { best = sc; bi = i; }
  }
  if (bi < 0) return p;
  p.a = bi + 3; p.e = e; p.extra = best;
  finalize(p, T[bidx(p.a, e)]);
  return p;
}

