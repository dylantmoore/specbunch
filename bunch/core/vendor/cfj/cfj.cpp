// cfj: engine for the CF.J fair test (notes/cfj_fair_test_scope_2026-11-09.md, "Operational
// definitions (pinned 2026-11-10)"). One pass per sample writes every band's fit and bootstrap
// p-values for
//   COS  on G0 (the incumbent LL.F3.S2 machinery), and
//   CF.J on G0, G1, G2,
// with p under both refit rules R0 (a non-converged bootstrap refit is dropped, as cfref) and
// R1 (it counts as an exceedance), from the same 499 multiplier draws. Window selection
// (S1 / S2 / S3) is glue over this table: py/cfj_select.py.
//
// Numerics are copied from the verified chain gqcsn.cpp -> llwin.cpp -> cfref.cpp
// (~/scratch-local/cjms-estimators/cppest/{spec,llwin,cftest/src}); the arithmetic of every
// fit, moment and bootstrap sum is unchanged. Changes, all of them bit-neutral relative to a
// single-threaded cfref run:
//   * Float reductions are fixed-order (no omp critical accumulation of doubles):
//       - cf_recurrence / sin_moments store per-512-block partial sums and add them in block
//         order (= the single-thread order of the originals);
//       - bootstrap multiplier sums are parallel over draws b, each draw's sum running over
//         the observations in order (= the single-thread order of boot_sums/boot_sums_sin).
//         The G0 columns use the originals' exact code shape (boot_g0); a reshaped version
//         changed the trig chain in the last bits, so shape matters here.
//     Results are therefore identical at any thread count.
//   * Grids have M points (G0 77, G1 97, G2 77..157); fixed-size scratch arrays are sized for
//     157. Moments and bootstrap sums are per grid column, and the longer grids reuse the G0
//     values for columns 0..76 (columns >= 77 come from boot_cols / the 157-point moment pass);
//     t_j and u_j use ugrid_ctx's formula with M points.
//   * Band tables are built in parallel over bands (each band's fit is independent).
//
// Grids (u = sigma_pilot * t, step .05 from .20): G0 [.20,4.00]; G1 [.20,5.00]; G2 [.20,u2],
// u2 = largest u on the .05 grid in [4.00,8.00] with m(t) >= SE(t), where m = sqrt((B E)^2 +
// (a E/t)^2), E = exp(-s^2 t^2/2), SE = sqrt((Vcc+Vss)/n), and (B, s, a) is the CF.J fit on the
// G0/R0/S2 rank-1 band ("the pilot"); u2 = 4.00 if none beyond 4.00 or no pilot.
// Joint information (for joint Info99 in the selector): per grid point, from the pilot under
// the 2x2 inverse-covariance (stacked) weights of fit_cfj's SE information,
//   IBB_j = w1 M_j^2,  ISS_j = w1 g^2 m1_j^2 + w2 g^2 m2_j^2,
//   g = s^2 t_j^2, M_j = B E_j, m1 = M_j - kap_j a x_j, m2 = a x_j, x_j = E_j / t_j
// (the (log B, log B) and (log s, log s) diagonal entries; gqcsn's Info99 sums the two).
//
// Usage: cfj LIST OUT.jsonl      LIST lines: path \t law \t n \t rep \t kink \t seed
// Output: one JSON line per sample (append mode), in list order:
//   {law, n, rep, seed, sigma_p,                      sigma_p = cosine Ensemble pilot (as cfref)
//    [internal: "refused:pilot" and nothing below but cpu, when the pilot fails]
//    pilot: {a, e, B, sigma, aJ} | null,              CF.J G0/R0/S2 rank-1 band and its fit
//    M2, u2,                                          G2 size and upper end
//    grids: {G0: {M, u[], t[], IBB[]|null, ISS[]|null, COS: TAB, CFJ: TAB},
//            G1: {M, u[], t[], IBB, ISS, CFJ: TAB},
//            G2: same as G1 | {same_as: "G0"} when u2 = 4.00},
//    cpu: {base, sin, boot_c, boot_s, u2, boot2, "COS.G0": [tab, boot_p], "CFJ.G0|G1|G2": [tab, boot_p], total},
//    wall, threads}
//   TAB = {n_valid, bands: [tuple, ...]}, n_valid = bands with a p-value. Band tuple
//   [a, e, B_raw, sigma, aJ, se_logB, se_logS, p_R0, p_R1, range_ok, converged, D]
//   a, e 0-based inclusive grid indices; every band of >= 11 points in (a, e) order; B/sigma/aJ/D
//   null unless converged; se null unless finite positive; p null unless converged with finite D
//   (cfref's candidate set); aJ null for COS; range_ok = llwin finalize rule. CPU values are
//   process CPU-s (all threads) of each stage. Numbers are %.17g (exact round trip).
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
#include <sys/resource.h>

using std::vector;
static const double NaN = std::numeric_limits<double>::quiet_NaN();
static double now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
static double cpu_now() {
  rusage r; getrusage(RUSAGE_SELF, &r);
  return r.ru_utime.tv_sec + 1e-6 * r.ru_utime.tv_usec + r.ru_stime.tv_sec + 1e-6 * r.ru_stime.tv_usec;
}
static std::string jnum(double x) {
  if (!std::isfinite(x)) return "null";
  char b[40]; std::snprintf(b, sizeof b, "%.17g", x); return b;
}

// ---------------------------------------------------------------- npy (gqcsn)
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

// ---------------------------------------------------------------- CF via recurrence (gqcsn; fixed-order)
// C[j] = mean cos(t_j dZ), C2[j] = mean cos^2, t_j = t0 + j*Delta, j=0..M-1.
static void cf_recurrence(const vector<double>& dZ, double t0, double Delta, int M,
                          vector<double>& C, vector<double>& C2) {
  const size_t n = dZ.size();
  C.assign(M, 0.0); C2.assign(M, 0.0);
  constexpr int BLK = 512;
  const size_t nblk = (n + BLK - 1) / BLK;
  vector<double> P1(nblk * (size_t)M), P2(nblk * (size_t)M);
#pragma omp parallel
  {
    alignas(64) double cc[BLK], ss[BLK], cd[BLK], sd[BLK];
#pragma omp for schedule(static)
    for (size_t b = 0; b < nblk; ++b) {
      size_t i0 = b * BLK, i1 = std::min(n, i0 + BLK); int L = (int)(i1 - i0);
      double* p1 = &P1[b * (size_t)M]; double* p2 = &P2[b * (size_t)M];
      for (int ii = 0; ii < L; ++ii) {
        double z = dZ[i0 + ii];
        cc[ii] = std::cos(t0 * z); ss[ii] = std::sin(t0 * z);
        cd[ii] = std::cos(Delta * z); sd[ii] = std::sin(Delta * z);
      }
      double a = 0, a2 = 0;
#pragma omp simd reduction(+:a,a2)
      for (int ii = 0; ii < L; ++ii) { a += cc[ii]; a2 += cc[ii] * cc[ii]; }
      p1[0] = a; p2[0] = a2;
      for (int j = 1; j < M; ++j) {
        a = 0; a2 = 0;
#pragma omp simd reduction(+:a,a2)
        for (int ii = 0; ii < L; ++ii) {
          double c = cc[ii], s = ss[ii];
          double cn = c * cd[ii] - s * sd[ii];
          double sn = s * cd[ii] + c * sd[ii];
          cc[ii] = cn; ss[ii] = sn; a += cn; a2 += cn * cn;
        }
        p1[j] = a; p2[j] = a2;
      }
    }
  }
  for (size_t b = 0; b < nblk; ++b)
    for (int j = 0; j < M; ++j) { C[j] += P1[b * (size_t)M + j]; C2[j] += P2[b * (size_t)M + j]; }
  const double inv_n = 1.0 / (double)n;
  for (int j = 0; j < M; ++j) { C[j] *= inv_n; C2[j] *= inv_n; }
}

// S = mean sin, S2 = mean sin^2, CS = mean cos*sin (cfref sin_moments; fixed-order).
static void sin_moments(const vector<double>& dZ, double t0, double Delta, int M,
                        vector<double>& S, vector<double>& S2, vector<double>& CS) {
  const size_t n = dZ.size();
  S.assign(M, 0.0); S2.assign(M, 0.0); CS.assign(M, 0.0);
  constexpr int BLK = 512;
  const size_t nblk = (n + BLK - 1) / BLK;
  vector<double> P1(nblk * (size_t)M), P2(nblk * (size_t)M), P3(nblk * (size_t)M);
#pragma omp parallel
  {
    alignas(64) double cc[BLK], ss[BLK], cd[BLK], sd[BLK];
#pragma omp for schedule(static)
    for (size_t b = 0; b < nblk; ++b) {
      size_t i0 = b * BLK, i1 = std::min(n, i0 + BLK); int L = (int)(i1 - i0);
      double* q1 = &P1[b * (size_t)M]; double* q2 = &P2[b * (size_t)M]; double* q3 = &P3[b * (size_t)M];
      for (int ii = 0; ii < L; ++ii) {
        double z = dZ[i0 + ii];
        cc[ii] = std::cos(t0 * z); ss[ii] = std::sin(t0 * z);
        cd[ii] = std::cos(Delta * z); sd[ii] = std::sin(Delta * z);
      }
      double s = 0, s2 = 0, cs = 0;
#pragma omp simd reduction(+:s,s2,cs)
      for (int ii = 0; ii < L; ++ii) { s += ss[ii]; s2 += ss[ii] * ss[ii]; cs += cc[ii] * ss[ii]; }
      q1[0] = s; q2[0] = s2; q3[0] = cs;
      for (int j = 1; j < M; ++j) {
        s = 0; s2 = 0; cs = 0;
#pragma omp simd reduction(+:s,s2,cs)
        for (int ii = 0; ii < L; ++ii) {
          double c = cc[ii], q = ss[ii];
          double cn = c * cd[ii] - q * sd[ii];
          double sn = q * cd[ii] + c * sd[ii];
          cc[ii] = cn; ss[ii] = sn; s += sn; s2 += sn * sn; cs += cn * sn;
        }
        q1[j] = s; q2[j] = s2; q3[j] = cs;
      }
    }
  }
  for (size_t b = 0; b < nblk; ++b)
    for (int j = 0; j < M; ++j) {
      S[j] += P1[b * (size_t)M + j]; S2[j] += P2[b * (size_t)M + j]; CS[j] += P3[b * (size_t)M + j];
    }
  const double inv_n = 1.0 / (double)n;
  for (int j = 0; j < M; ++j) { S[j] *= inv_n; S2[j] *= inv_n; CS[j] *= inv_n; }
}

// ---------------------------------------------------------------- Ensemble sigma pilot (gqcsn, verbatim)
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

// ---------------------------------------------------------------- grids
constexpr double UGRID_LO = 0.20, UGRID_DU = 0.05;
constexpr int NG0 = 77, NG1 = 97, NGMAX = 157, NGPAD = 160;  // u <= 4.00, 5.00, 8.00
constexpr int MP0 = 80;     // phase-1 bootstrap columns (G0's padded row, as gqcsn MP)
constexpr int BOOT_B = 499;
constexpr int MINPTS = 11;

struct Ctx {
  size_t n; int M; double sigma_p, t_start, Delta;
  vector<double> u, t, C, var, s;
};
// ugrid_ctx with M points, from moments already computed on the 157-point recurrence grid.
static Ctx make_ctx(size_t n, double sp, int M, const vector<double>& Cf, const vector<double>& C2f) {
  Ctx c; c.n = n; c.M = M; c.sigma_p = sp;
  c.t_start = UGRID_LO / sp; c.Delta = UGRID_DU / sp;
  const double UHI = (20.0 + 5.0 * (M - 1)) / 100.0;
  c.C.assign(Cf.begin(), Cf.begin() + M);
  c.var.resize(M); c.t.resize(M); c.u.resize(M); c.s.resize(M);
  double tend = c.t_start + (M - 1) * c.Delta;
  for (int j = 0; j < M; ++j) {
    c.var[j] = C2f[j] - c.C[j] * c.C[j];
    c.t[j] = (j == M - 1) ? tend : c.t_start + j * ((tend - c.t_start) / (M - 1));
    c.u[j] = (j == M - 1) ? UHI : UGRID_LO + j * ((UHI - UGRID_LO) / (M - 1));
    c.s[j] = std::sqrt(std::max(c.var[j] / (double)c.n, 1e-20));
  }
  return c;
}

// ---------------------------------------------------------------- fits
struct BF {
  double B = NaN, sigma = NaN, logB = NaN, logS = NaN, seB = NaN, seS = NaN, D = NaN;
  bool valid = false, sevalid = false;
};

// gqcsn gqml_nr (verbatim but for the scratch size)
struct NR { double B, sigma; bool conv; };
static NR gqml_nr(const double* C, const double* t, const double* var, int M, size_t n,
                  double sigma_init, int max_iter = 100, double tol = 1e-8) {
  double w[NGPAD];
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

// llwin fit_lev (F3): gqml_nr, SE from Gauss-Newton weighted information in (log B, log sigma).
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

// Per-grid sine moments and CF.J weights (cfref make_mom, the COS/CF parts).
struct Mom { vector<double> S, Vss, Vcs, kap, varCt, sS, w1, w2; };
static Mom make_mom(const Ctx& c, const vector<double>& S, const vector<double>& S2, const vector<double>& CS) {
  Mom M; const double nd = (double)c.n; const int NG = c.M;
  M.S.assign(S.begin(), S.begin() + NG); M.Vss.resize(NG); M.Vcs.resize(NG);
  M.kap.resize(NG); M.varCt.resize(NG); M.sS.resize(NG); M.w1.resize(NG); M.w2.resize(NG);
  for (int j = 0; j < NG; ++j) {
    const double C = c.C[j], Sj = S[j], Vcc = c.var[j];
    M.Vss[j] = S2[j] - Sj * Sj; M.Vcs[j] = CS[j] - C * Sj;
    M.kap[j] = M.Vcs[j] / M.Vss[j];
    M.varCt[j] = Vcc - M.Vcs[j] * M.Vcs[j] / M.Vss[j];
    M.sS[j] = std::sqrt(std::max(M.Vss[j] / nd, 1e-20));
    M.w1[j] = M.varCt[j] > 1e-12 ? nd / M.varCt[j] : nd;
    M.w2[j] = M.Vss[j] > 1e-12 ? nd / M.Vss[j] : nd;
  }
  return M;
}

// cfref fit_cfj (verbatim but for the scratch size)
template <bool WANT_SE>
static bool fit_cfj(const double* Cv, const double* Sv, const double* t, const double* kap, const double* w1,
                    const double* w2, int m, double sigma_init, BF& r, double& aout) {
  double Ct[NGPAD];
  for (int j = 0; j < m; ++j) Ct[j] = Cv[j] - kap[j] * Sv[j];
  double sxx = 0, sx = 0, sy = 0, sxy = 0;
  for (int j = 0; j < m; ++j) {
    double x = -t[j] * t[j] / 2.0, y = std::log(std::max(Cv[j], 1e-10));
    sx += x; sy += y; sxx += x * x; sxy += x * y;
  }
  double Md = (double)m, denom = Md * sxx - sx * sx, al, s2i;
  if (std::fabs(denom) > 1e-30) {
    double b2 = (Md * sxy - sx * sy) / denom;
    al = (sy - b2 * sx) / Md;
    s2i = std::max((Md * sxy - sx * sy) / denom, sigma_init * sigma_init);
  } else { al = std::log(std::max(Cv[0], 1e-10)); s2i = sigma_init * sigma_init; }
  double ps = 0.5 * std::log(s2i);
  double aa;
  {
    double s2 = std::exp(2.0 * ps), num = 0, den = 0;
    for (int j = 0; j < m; ++j) {
      double g = s2 * t[j] * t[j], E = std::exp(std::clamp(-g / 2.0, -500.0, 500.0)), x = E / t[j];
      double Mj = std::exp(std::clamp(al - g / 2.0, -500.0, 500.0));
      num += w2[j] * Sv[j] * x - w1[j] * kap[j] * x * (Ct[j] - Mj);
      den += (w2[j] + w1[j] * kap[j] * kap[j]) * x * x;
    }
    aa = den > 0 ? num / den : 0.0;
  }
  bool conv = false;
  for (int it = 0; it < 100; ++it) {
    double s2 = std::exp(2.0 * ps);
    double gr[3] = {0, 0, 0}, H[3][3] = {{0}}, G[3][3] = {{0}};
    for (int j = 0; j < m; ++j) {
      double g = s2 * t[j] * t[j], E = std::exp(std::clamp(-g / 2.0, -500.0, 500.0)), x = E / t[j];
      double Mj = std::exp(std::clamp(al - g / 2.0, -500.0, 500.0));
      double m1 = Mj - kap[j] * aa * x, r1 = Ct[j] - m1;
      double m2 = aa * x, r2 = Sv[j] - m2;
      double d1[3] = {Mj, -g * m1, -kap[j] * x}, d2[3] = {0.0, -g * m2, x};
      double h1[3][3] = {{Mj, -g * Mj, 0.0}, {-g * Mj, (g * g - 2.0 * g) * m1, kap[j] * g * x}, {0.0, kap[j] * g * x, 0.0}};
      double h2[3][3] = {{0.0, 0.0, 0.0}, {0.0, (g * g - 2.0 * g) * m2, -g * x}, {0.0, -g * x, 0.0}};
      for (int p = 0; p < 3; ++p) {
        gr[p] += w1[j] * r1 * d1[p] + w2[j] * r2 * d2[p];
        for (int q = 0; q < 3; ++q) {
          double gg = -(w1[j] * d1[p] * d1[q] + w2[j] * d2[p] * d2[q]);
          G[p][q] += gg; H[p][q] += gg + w1[j] * r1 * h1[p][q] + w2[j] * r2 * h2[p][q];
        }
      }
    }
    auto det2 = [](double A[3][3]) { return A[0][0] * A[1][1] - A[0][1] * A[1][0]; };
    auto det3 = [](double A[3][3]) {
      return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
             A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]); };
    double (*A)[3];
    if (H[0][0] < 0.0 && det2(H) > 1e-30 && det3(H) < -1e-30) A = H; else A = G;
    double d = det3(A);
    if (std::fabs(d) < 1e-30) break;
    double inv[3][3];
    inv[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / d; inv[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) / d;
    inv[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / d; inv[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) / d;
    inv[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / d; inv[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) / d;
    inv[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) / d; inv[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1]) / d;
    inv[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / d;
    double dl[3];
    for (int p = 0; p < 3; ++p) dl[p] = inv[p][0] * gr[0] + inv[p][1] * gr[1] + inv[p][2] * gr[2];
    if (!std::isfinite(dl[0]) || !std::isfinite(dl[1]) || !std::isfinite(dl[2])) break;
    al -= dl[0]; ps -= dl[1]; aa -= dl[2];
    if (dl[0] * dl[0] + dl[1] * dl[1] + dl[2] * dl[2] < 1e-16 * (1.0 + al * al + ps * ps + aa * aa)) { conv = true; break; }
  }
  const double B = std::exp(std::clamp(al, -500.0, 500.0)), sg = std::exp(ps);
  if (!(conv && std::isfinite(B) && std::isfinite(sg) && std::isfinite(aa))) return false;
  r.B = B; r.sigma = sg; r.logB = std::log(B); r.logS = std::log(sg); aout = aa;
  if (WANT_SE) {
    const double s2 = sg * sg;
    double I[3][3] = {{0}};
    for (int j = 0; j < m; ++j) {
      double g = s2 * t[j] * t[j], E = std::exp(-g / 2.0), x = E / t[j], Mj = B * E;
      double m1 = Mj - kap[j] * aa * x, m2 = aa * x;
      double d1[3] = {Mj, -g * m1, -kap[j] * x}, d2[3] = {0.0, -g * m2, x};
      for (int p = 0; p < 3; ++p) for (int q = 0; q < 3; ++q) I[p][q] += w1[j] * d1[p] * d1[q] + w2[j] * d2[p] * d2[q];
    }
    double d = I[0][0] * (I[1][1] * I[2][2] - I[1][2] * I[2][1]) - I[0][1] * (I[1][0] * I[2][2] - I[1][2] * I[2][0]) +
               I[0][2] * (I[1][0] * I[2][1] - I[1][1] * I[2][0]);
    double v00 = (I[1][1] * I[2][2] - I[1][2] * I[2][1]) / d, v11 = (I[0][0] * I[2][2] - I[0][2] * I[2][0]) / d;
    r.seB = std::sqrt(v00); r.seS = std::sqrt(v11);
    r.sevalid = d > 0 && std::isfinite(r.seB) && std::isfinite(r.seS) && r.seB > 0 && r.seS > 0;
  }
  return true;
}

// Fit + sup standardized residual (cfref vfit, variants COS = 0 and CF.J = 4 only).
struct VF { BF r; double p1 = NaN; };
template <bool WANT_SE>
static bool vfit(bool cfj, const Ctx& c, const Mom& M, const double* Y1, const double* Y2, int a, int m,
                 double sinit, VF& o) {
  const double* t = c.t.data() + a;
  if (!cfj) {
    bool ok = fit_lev<WANT_SE>(Y1 + a, c.t.data() + a, c.var.data() + a, m, c.n, sinit, o.r);
    if (ok) o.r.D = sup_std_resid(&Y1[a], &c.t[a], &c.s[a], m, o.r.B, o.r.sigma);
    return ok;
  }
  double aj = NaN;
  bool ok = fit_cfj<WANT_SE>(Y1 + a, Y2 + a, t, M.kap.data() + a, M.w1.data() + a, M.w2.data() + a, m, sinit, o.r, aj);
  if (!ok) return false;
  o.p1 = aj;
  double D = 0;
  const double s2 = o.r.sigma * o.r.sigma;
  for (int j = 0; j < m; ++j) {
    double E = std::exp(-s2 * t[j] * t[j] * 0.5);
    D = std::max(D, std::fabs(Y1[a + j] - o.r.B * E) / c.s[a + j]);
    D = std::max(D, std::fabs(Y2[a + j] - aj * E / t[j]) / M.sS[a + j]);
  }
  o.r.D = D;
  return true;
}

// ---------------------------------------------------------------- bootstrap multiplier sums
static inline uint64_t splitmix64(uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}
static inline uint64_t sign_word(uint64_t key, int b, uint64_t g) {
  return splitmix64(key ^ splitmix64(((uint64_t)b << 40) ^ g));
}

// G0 columns (phase 1): gqcsn boot_sums (SIN = false; also S[b] = sum_i xi_{b,i}) and cfref
// boot_sums_sin (SIN = true), code shape kept exactly (X[64][80], constant 77-point chain): the
// compiler's evaluation of the trig chain depends on that shape, and these sums must match
// cfref bitwise. Only the parallelization changed: over draws b instead of observation groups
// with an omp critical sum, so each (b, j) sum runs over the observations in order (= the
// single-threaded original) at any thread count. Y has row length MP0 = 80.
constexpr int UGRID_N = NG0, MP = MP0;
template <bool SIN>
static void boot_g0(const vector<double>& dZ, double t0, double Delta, uint64_t seed,
                    vector<double>& Y, vector<double>* S) {
  const size_t n = dZ.size();
  const uint64_t key = splitmix64((uint64_t)seed + 7777777ULL);
  const size_t ng = (n + 63) / 64;
  Y.assign((size_t)BOOT_B * MP, 0.0);
  if (S) S->assign(BOOT_B, 0.0);
#pragma omp parallel
  {
    const int nth = omp_get_num_threads(), tid = omp_get_thread_num();
    const int b0 = (int)((long)BOOT_B * tid / nth), b1 = (int)((long)BOOT_B * (tid + 1) / nth);
    alignas(64) double X[64][MP];
    if (b1 > b0)
    for (size_t g = 0; g < ng; ++g) {
      size_t i0 = g * 64; int L = (int)std::min<size_t>(64, n - i0);
      for (int ii = 0; ii < 64; ++ii) {
        if (ii >= L) { for (int j = 0; j < MP; ++j) X[ii][j] = 0.0; continue; }
        double z = dZ[i0 + ii];
        double cl = std::cos(t0 * z), sl = std::sin(t0 * z), cd = std::cos(Delta * z), sd = std::sin(Delta * z);
        if (!SIN) {
          X[ii][0] = cl; double cp = cl, sp = sl;
          for (int j = 1; j < UGRID_N; ++j) {
            double cn = cp * cd - sp * sd, sn = sp * cd + cp * sd;
            X[ii][j] = cn; cp = cn; sp = sn;
          }
        } else {
          X[ii][0] = sl; double cp = cl, sp = sl;
          for (int j = 1; j < UGRID_N; ++j) {
            double cn = cp * cd - sp * sd, sn = sp * cd + cp * sd;
            X[ii][j] = sn; cp = cn; sp = sn;
          }
        }
        for (int j = UGRID_N; j < MP; ++j) X[ii][j] = 0.0;
      }
      const uint64_t mask = (L == 64) ? ~0ULL : ((1ULL << L) - 1);
      for (int b = b0; b < b1; ++b) {
        uint64_t w = sign_word(key, b, g);
        if (S) (*S)[b] += 2.0 * __builtin_popcountll(w & mask) - (double)L;
        double* yb = &Y[(size_t)b * MP];
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

// Extension columns (phase 2, G1/G2): columns [j0, j1) (multiples of 16, j1 <= NGPAD) of Y[b*NGPAD + j] = sum_i xi_{b,i} f(t_j dZ_i),
// f = cos (SIN = false, gqcsn boot_sums) or sin (cfref boot_sums_sin); S[b] = sum_i xi_{b,i} when
// S != nullptr. Same multipliers and, per (b, j), the same observation-order accumulation as the
// single-threaded originals. Parallel over draws; each thread recomputes the trig rows.
template <bool SIN>
static void boot_cols(const vector<double>& dZ, double t0, double Delta, uint64_t seed, int j0, int j1,
                      vector<double>& Y, vector<double>* S) {
  const size_t n = dZ.size();
  const uint64_t key = splitmix64((uint64_t)seed + 7777777ULL);
  const size_t ng = (n + 63) / 64;
  const int jr = std::min(j1, NGMAX);  // real grid columns
  if (Y.size() != (size_t)BOOT_B * NGPAD) Y.assign((size_t)BOOT_B * NGPAD, 0.0);
  if (S) S->assign(BOOT_B, 0.0);
#pragma omp parallel
  {
    const int nth = omp_get_num_threads(), tid = omp_get_thread_num();
    const int b0 = (int)((long)BOOT_B * tid / nth), b1 = (int)((long)BOOT_B * (tid + 1) / nth);
    vector<double> Xv((size_t)64 * NGPAD, 0.0);
    double* X = Xv.data();
    if (b1 > b0)
      for (size_t g = 0; g < ng; ++g) {
        size_t i0 = g * 64; int L = (int)std::min<size_t>(64, n - i0);
        for (int ii = 0; ii < 64; ++ii) {
          double* Xi = X + (size_t)ii * NGPAD;
          if (ii >= L) { for (int j = 0; j < j1; ++j) Xi[j] = 0.0; continue; }
          double z = dZ[i0 + ii];
          double cl = std::cos(t0 * z), sl = std::sin(t0 * z), cd = std::cos(Delta * z), sd = std::sin(Delta * z);
          Xi[0] = SIN ? sl : cl; double cp = cl, sp = sl;
          for (int j = 1; j < jr; ++j) {
            double cn = cp * cd - sp * sd, sn = sp * cd + cp * sd;
            Xi[j] = SIN ? sn : cn; cp = cn; sp = sn;
          }
          for (int j = jr; j < j1; ++j) Xi[j] = 0.0;
        }
        const uint64_t mask = (L == 64) ? ~0ULL : ((1ULL << L) - 1);
        for (int b = b0; b < b1; ++b) {
          uint64_t w = sign_word(key, b, g);
          if (S) (*S)[b] += 2.0 * __builtin_popcountll(w & mask) - (double)L;
          double* yb = &Y[(size_t)b * NGPAD];
          for (int jc = j0; jc < j1; jc += 16) {
            double acc[16];
            for (int q = 0; q < 16; ++q) acc[q] = yb[jc + q];
            for (int ii = 0; ii < 64; ++ii) {
              double s = ((w >> ii) & 1ULL) ? 1.0 : -1.0;
              const double* Xi = X + (size_t)ii * NGPAD;
#pragma omp simd
              for (int q = 0; q < 16; ++q) acc[q] += s * Xi[jc + q];
            }
            for (int q = 0; q < 16; ++q) yb[jc + q] = acc[q];
          }
        }
      }
  }
}

// ---------------------------------------------------------------- one (grid, variant) pass
struct Table {
  int M = 0; bool cfj = false;
  vector<VF> T; vector<double> pR0, pR1; vector<int> ids;  // ids: valid bands with finite D
  double t_tab = 0, t_s2 = 0;
};

static Table run_variant(bool cfj, const Ctx& c, const Mom& M, const vector<double>& Yc, const vector<double>& Ys,
                         const vector<double>& Sb) {
  Table R; R.M = c.M; R.cfj = cfj;
  const int NG = c.M;
  auto bidx = [NG](int a, int e) { return a * NG + e; };
  const double* Y1 = c.C.data();
  const double* Y2 = M.S.data();
  const double sinit = c.sigma_p;
  double c0 = cpu_now();
  R.T.assign((size_t)NG * NG, VF());
  vector<double> MU((size_t)NG * NG * NG, 0.0), MU2;
  if (cfj) MU2.assign((size_t)NG * NG * NG, 0.0);
  vector<int> all;
  for (int a = 0; a < NG; ++a)
    for (int e = a + MINPTS - 1; e < NG; ++e) all.push_back(bidx(a, e));
#pragma omp parallel for schedule(dynamic, 8)
  for (size_t q = 0; q < all.size(); ++q) {
    const int id = all[q], a = id / NG, e = id % NG, m = e - a + 1;
    VF o;
    bool ok = vfit<true>(cfj, c, M, Y1, Y2, a, m, sinit, o);
    if (ok) {
      double* mu = &MU[(size_t)id * NG];
      for (int j = a; j <= e; ++j) mu[j] = mu_k0(o.r.B, o.r.sigma, c.t[j]);
      if (cfj) {
        double* mu2 = &MU2[(size_t)id * NG];
        for (int j = a; j <= e; ++j) mu2[j] = o.p1 * std::exp(-o.r.sigma * o.r.sigma * c.t[j] * c.t[j] * 0.5) / c.t[j];
      }
    }
    o.r.valid = ok;
    R.T[id] = o;
  }
  R.t_tab = cpu_now() - c0;
  c0 = cpu_now();
  for (int id : all)
    if (R.T[id].r.valid && std::isfinite(R.T[id].r.D)) R.ids.push_back(id);
  const int nb = (int)R.ids.size();
  // per-draw exceedance / failure flags, summed in draw order afterwards (integers: exact anyway)
  vector<int> exc((size_t)nb, 0), fail((size_t)nb, 0);
  const double inv_n = 1.0 / (double)c.n;
#pragma omp parallel
  {
    vector<int> el((size_t)nb, 0), fl((size_t)nb, 0);
    double Cs[NGPAD], Ss[NGPAD], dc[NGPAD], ds[NGPAD];
#pragma omp for schedule(dynamic, 4)
    for (int b = 0; b < BOOT_B; ++b) {
      const double* yb = &Yc[(size_t)b * NGPAD];
      const double* ysb = &Ys[(size_t)b * NGPAD];
      const double Sbb = Sb[b];
      if (cfj)
        for (int j = 0; j < NG; ++j) { dc[j] = (yb[j] - c.C[j] * Sbb) * inv_n; ds[j] = (ysb[j] - M.S[j] * Sbb) * inv_n; }
      for (int k = 0; k < nb; ++k) {
        const int id = R.ids[k], a = id / NG, e = id % NG, m = e - a + 1;
        const double* mub = &MU[(size_t)id * NG];
        const VF& tb = R.T[id];
        if (!cfj) {
          for (int j = a; j <= e; ++j) Cs[j] = mub[j] + (yb[j] - c.C[j] * Sbb) * inv_n;  // llref expression
        } else {
          const double* mu2 = &MU2[(size_t)id * NG];
          for (int j = a; j <= e; ++j) { Cs[j] = mub[j] + dc[j]; Ss[j] = mu2[j] + ds[j]; }
        }
        VF o;
        if (!vfit<false>(cfj, c, M, Cs, Ss, a, m, tb.r.sigma, o)) { fl[k]++; continue; }
        if (o.r.D >= tb.r.D) el[k]++;
      }
    }
#pragma omp critical
    for (int k = 0; k < nb; ++k) { exc[k] += el[k]; fail[k] += fl[k]; }
  }
  R.pR0.assign((size_t)NG * NG, NaN); R.pR1.assign((size_t)NG * NG, NaN);
  for (int k = 0; k < nb; ++k) {
    const int id = R.ids[k];
    R.pR0[id] = (1.0 + exc[k]) / (BOOT_B + 1.0);
    R.pR1[id] = (1.0 + (exc[k] + fail[k])) / (BOOT_B + 1.0);
  }
  R.t_s2 = cpu_now() - c0;
  return R;
}

// cfref S2 ranking (p >= .10 by most points, ties lower u_L; then the rest by descending p, ties
// most points, lower u_L). Used here only to find the CF.J G0/R0/S2 pilot band.
static vector<int> s2_rank(const Table& R, const vector<double>& PV) {
  const int NG = R.M;
  vector<int> hi, lo;
  for (int id : R.ids) (PV[id] >= 0.10 ? hi : lo).push_back(id);
  auto pts = [NG](int id) { return id % NG - id / NG + 1; };
  std::stable_sort(hi.begin(), hi.end(), [&](int x, int y) {
    if (pts(x) != pts(y)) return pts(x) > pts(y);
    return x / NG < y / NG; });
  std::stable_sort(lo.begin(), lo.end(), [&](int x, int y) {
    if (PV[x] != PV[y]) return PV[x] > PV[y];
    if (pts(x) != pts(y)) return pts(x) > pts(y);
    return x / NG < y / NG; });
  vector<int> rk(hi); rk.insert(rk.end(), lo.begin(), lo.end());
  return rk;
}

static std::string arr(const vector<double>& x) {
  std::string s = "[";
  for (size_t i = 0; i < x.size(); ++i) s += (i ? "," : "") + jnum(x[i]);
  return s + "]";
}

static std::string bands_json(const Table& R) {
  const int NG = R.M;
  std::string s = "{\"n_valid\":" + std::to_string(R.ids.size()) + ",\"bands\":[";
  bool first = true;
  for (int a = 0; a < NG; ++a)
    for (int e = a + MINPTS - 1; e < NG; ++e) {
      const int id = a * NG + e; const VF& o = R.T[id]; const BF& r = o.r;
      bool rok = r.valid && std::isfinite(r.B) && r.B > 0.0 && r.B <= 1.5 && std::isfinite(r.sigma) && r.sigma >= 1e-3 && r.sigma <= 5.0;
      s += std::string(first ? "" : ",") + "[" + std::to_string(a) + "," + std::to_string(e) + "," +
           (r.valid ? jnum(r.B) : "null") + "," + (r.valid ? jnum(r.sigma) : "null") + "," +
           (r.valid && R.cfj ? jnum(o.p1) : "null") + "," +
           (r.valid && r.sevalid ? jnum(r.seB) : "null") + "," + (r.valid && r.sevalid ? jnum(r.seS) : "null") + "," +
           jnum(R.pR0[id]) + "," + jnum(R.pR1[id]) + "," + (rok ? "1" : "0") + "," + (r.valid ? "1" : "0") + "," +
           (r.valid ? jnum(r.D) : "null") + "]";
      first = false;
    }
  return s + "]}";
}

// Joint per-point information from the pilot (B, s, a) on grid c (see header).
static void joint_info(const Ctx& c, const Mom& M, double B, double sg, double aa, vector<double>& IBB, vector<double>& ISS) {
  const int NG = c.M; IBB.assign(NG, NaN); ISS.assign(NG, NaN);
  const double s2 = sg * sg;
  for (int j = 0; j < NG; ++j) {
    const double t = c.t[j], g = s2 * t * t, E = std::exp(-g / 2.0), x = E / t, Mj = B * E;
    const double m1 = Mj - M.kap[j] * aa * x, m2 = aa * x;
    IBB[j] = M.w1[j] * Mj * Mj;
    ISS[j] = M.w1[j] * (g * m1) * (g * m1) + M.w2[j] * (g * m2) * (g * m2);
  }
}

int main(int argc, char** argv) {
  if (argc != 3) { std::fprintf(stderr, "usage: cfj LIST OUT.jsonl\n"); return 2; }
  std::ifstream lf(argv[1]);
  if (!lf) { std::perror(argv[1]); return 1; }
  std::FILE* out = std::fopen(argv[2], "a");
  if (!out) { std::perror(argv[2]); return 1; }
  std::string line;
  while (std::getline(lf, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line); std::string path, law, ns, reps, ks, seeds;
    std::getline(ss, path, '\t'); std::getline(ss, law, '\t'); std::getline(ss, ns, '\t');
    std::getline(ss, reps, '\t'); std::getline(ss, ks, '\t'); std::getline(ss, seeds, '\t');
    double k = std::atof(ks.c_str()); uint64_t seed = std::strtoull(seeds.c_str(), nullptr, 10);
    vector<double> Z;
    try { Z = read_npy(path); } catch (std::exception& ex) { std::fprintf(stderr, "READFAIL %s %s\n", path.c_str(), ex.what()); continue; }
    const size_t n = Z.size();
    const double W0 = now(), P0 = cpu_now();
    vector<double> dZ(n);
    for (size_t i = 0; i < n; ++i) dZ[i] = Z[i] - k;
    vector<double>().swap(Z);
    std::string dj = "{\"law\":\"" + law + "\",\"n\":" + std::to_string(n) + ",\"rep\":" + reps + ",\"seed\":" + seeds;
    Pilot P = ensemble_sigma(dZ);
    dj += ",\"sigma_p\":" + jnum(P.sigma_p);
    double t_base = NaN, t_sin = 0, t_bc = 0, t_bs = 0, t_b2 = 0, t_u2 = 0;
    Table TC0, TJ[3];
    std::string cpuv;
    if (!(std::isfinite(P.sigma_p) && P.sigma_p > 0)) {
      t_base = cpu_now() - P0;
      dj += ",\"internal\":\"refused:pilot\"";
    } else {
      const double t_start = UGRID_LO / P.sigma_p, Delta = UGRID_DU / P.sigma_p;
      vector<double> Cf, C2f;
      cf_recurrence(dZ, t_start, Delta, NGMAX, Cf, C2f);
      Ctx cg[3];
      cg[0] = make_ctx(n, P.sigma_p, NG0, Cf, C2f);
      t_base = cpu_now() - P0;
      double c0 = cpu_now();
      vector<double> S, S2, CS;
      sin_moments(dZ, t_start, Delta, NGMAX, S, S2, CS);
      Mom Mg[3];
      Mg[0] = make_mom(cg[0], S, S2, CS);
      t_sin = cpu_now() - c0;
      // phase-1 bootstrap sums (columns 0..79), cos then sin; t0/Delta exactly as cfref passes them
      vector<double> Yc0, Ys0, Sb, Yc((size_t)BOOT_B * NGPAD, 0.0), Ys((size_t)BOOT_B * NGPAD, 0.0);
      c0 = cpu_now();
      boot_g0<false>(dZ, cg[0].t[0], cg[0].t[1] - cg[0].t[0], seed, Yc0, &Sb);
      for (int b = 0; b < BOOT_B; ++b) std::copy(&Yc0[(size_t)b * MP0], &Yc0[(size_t)b * MP0] + NG0, &Yc[(size_t)b * NGPAD]);
      t_bc = cpu_now() - c0;
      c0 = cpu_now();
      boot_g0<true>(dZ, cg[0].t[0], cg[0].t[1] - cg[0].t[0], seed, Ys0, nullptr);
      for (int b = 0; b < BOOT_B; ++b) std::copy(&Ys0[(size_t)b * MP0], &Ys0[(size_t)b * MP0] + NG0, &Ys[(size_t)b * NGPAD]);
      t_bs = cpu_now() - c0;
      // G0: COS and CF.J
      TC0 = run_variant(false, cg[0], Mg[0], Yc, Ys, Sb);
      TJ[0] = run_variant(true, cg[0], Mg[0], Yc, Ys, Sb);
      // pilot: CF.J G0 / R0 / S2 rank-1 band
      c0 = cpu_now();
      vector<int> rk = s2_rank(TJ[0], TJ[0].pR0);
      bool havep = !rk.empty();
      double pB = NaN, pS = NaN, pA = NaN; int pa = -1, pe = -1;
      if (havep) { const VF& o = TJ[0].T[rk[0]]; pB = o.r.B; pS = o.r.sigma; pA = o.p1; pa = rk[0] / NG0; pe = rk[0] % NG0; }
      // u2 on the 157-point grid: largest j in [76, 156] with m >= SE; else 76
      Ctx cmax = make_ctx(n, P.sigma_p, NGMAX, Cf, C2f);
      int j2 = NG0 - 1;
      if (havep) {
        for (int j = NG0 - 1; j < NGMAX; ++j) {
          const double t = cmax.t[j], E = std::exp(-pS * pS * t * t / 2.0);
          const double mB = pB * E, mA = pA * E / t;
          const double m = std::sqrt(mB * mB + mA * mA);
          const double Vss = S2[j] - S[j] * S[j];
          const double se = std::sqrt((cmax.var[j] + Vss) / (double)n);
          if (m >= se) j2 = j;
        }
      }
      const int M2 = j2 + 1;
      t_u2 = cpu_now() - c0;
      dj += ",\"pilot\":" + (havep ? "{\"a\":" + std::to_string(pa) + ",\"e\":" + std::to_string(pe) + ",\"B\":" + jnum(pB) +
                                       ",\"sigma\":" + jnum(pS) + ",\"aJ\":" + jnum(pA) + "}" : std::string("null"));
      dj += ",\"M2\":" + std::to_string(M2) + ",\"u2\":" + jnum((20.0 + 5.0 * (M2 - 1)) / 100.0);
      // phase-2 bootstrap columns: computed for [64, W), columns >= 77 copied in (0..76 stay G0's)
      const int W = std::min(NGPAD, ((std::max(NG1, M2) + 15) / 16) * 16);
      c0 = cpu_now();
      {
        vector<double> Ye;
        boot_cols<false>(dZ, cg[0].t[0], cg[0].t[1] - cg[0].t[0], seed, 64, W, Ye, nullptr);
        for (int b = 0; b < BOOT_B; ++b) std::copy(&Ye[(size_t)b * NGPAD + NG0], &Ye[(size_t)b * NGPAD + W], &Yc[(size_t)b * NGPAD + NG0]);
        Ye.clear();
        boot_cols<true>(dZ, cg[0].t[0], cg[0].t[1] - cg[0].t[0], seed, 64, W, Ye, nullptr);
        for (int b = 0; b < BOOT_B; ++b) std::copy(&Ye[(size_t)b * NGPAD + NG0], &Ye[(size_t)b * NGPAD + W], &Ys[(size_t)b * NGPAD + NG0]);
      }
      t_b2 = cpu_now() - c0;
      cg[1] = make_ctx(n, P.sigma_p, NG1, Cf, C2f); Mg[1] = make_mom(cg[1], S, S2, CS);
      TJ[1] = run_variant(true, cg[1], Mg[1], Yc, Ys, Sb);
      const bool g2same = (M2 == NG0);
      if (!g2same) {
        cg[2] = make_ctx(n, P.sigma_p, M2, Cf, C2f); Mg[2] = make_mom(cg[2], S, S2, CS);
        TJ[2] = run_variant(true, cg[2], Mg[2], Yc, Ys, Sb);
      }
      // output
      dj += ",\"grids\":{";
      for (int gi = 0; gi < 3; ++gi) {
        dj += std::string(gi ? "," : "") + "\"G" + std::to_string(gi) + "\":";
        if (gi == 2 && g2same) { dj += "{\"same_as\":\"G0\"}"; continue; }
        const Ctx& c = cg[gi];
        vector<double> IBB, ISS;
        if (havep) joint_info(c, Mg[gi], pB, pS, pA, IBB, ISS);
        dj += "{\"M\":" + std::to_string(c.M) + ",\"u\":" + arr(c.u) + ",\"t\":" + arr(c.t) +
              ",\"IBB\":" + (havep ? arr(IBB) : std::string("null")) + ",\"ISS\":" + (havep ? arr(ISS) : std::string("null"));
        if (gi == 0) dj += ",\"COS\":" + bands_json(TC0);
        dj += ",\"CFJ\":" + bands_json(TJ[gi]) + "}";
      }
      dj += "}";
      cpuv = ",\"COS.G0\":[" + jnum(TC0.t_tab) + "," + jnum(TC0.t_s2) + "]";
      for (int gi = 0; gi < 3; ++gi)
        cpuv += ",\"CFJ.G" + std::to_string(gi) + "\":[" + jnum(gi == 2 && g2same ? TJ[0].t_tab : TJ[gi].t_tab) + "," +
                jnum(gi == 2 && g2same ? TJ[0].t_s2 : TJ[gi].t_s2) + "]";
    }
    double cpu = cpu_now() - P0;
    dj += ",\"cpu\":{\"base\":" + jnum(t_base) + ",\"sin\":" + jnum(t_sin) + ",\"boot_c\":" + jnum(t_bc) + ",\"boot_s\":" + jnum(t_bs) +
          ",\"u2\":" + jnum(t_u2) + ",\"boot2\":" + jnum(t_b2) + cpuv + ",\"total\":" + jnum(cpu) + "},\"wall\":" + jnum(now() - W0) +
          ",\"threads\":" + std::to_string(omp_get_max_threads()) + "}\n";
    std::fputs(dj.c_str(), out); std::fflush(out);
    std::fprintf(stderr, "%s n=%zu rep=%s wall=%.2fs cpu=%.2fs\n", law.c_str(), n, reps.c_str(), now() - W0, cpu);
  }
  std::fclose(out);
  return 0;
}
