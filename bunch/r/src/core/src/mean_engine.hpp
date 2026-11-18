




































#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <utility>
#include <vector>
#include "p2/p2pilot.hpp"

namespace p2 {
using std::vector;
static const double NaN = std::numeric_limits<double>::quiet_NaN();
constexpr double PI = 3.141592653589793238462643383279502884;

static inline uint64_t sm64(uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}
struct CStream {  
  uint64_t base, k = 0;
  CStream(uint64_t seed, uint64_t i) : base(sm64(seed ^ sm64(i + 0x632BE59BD9B4E019ULL))) {}
  double u() { return ((double)(sm64(base + (k++)) >> 12) + 0.5) * 0x1p-52; }
};
static double Qn(double u) {  
  const double SQ2 = 1.414213562373095048801688724209698079;
  double p = std::min(u, 1 - u), x;
  if (p < 0.02425) {
    double t = std::sqrt(-2 * std::log(p));
    x = -(((((-.007784894002430293 * t - .3223964580411365) * t - 2.400758277161838) * t - 2.549732539343734) * t + 4.374664141464968) * t + 2.938163982698783) /
        ((((.007784695709041462 * t + .3224671290700398) * t + 2.445134137142996) * t + 3.754408661907416) * t + 1);
  } else {
    double q = p - .5, r = q * q;
    x = -(((((-39.69683028665376 * r + 220.9460984245205) * r - 275.9285104469687) * r + 138.3577518672690) * r - 30.66479806614716) * r + 2.506628277459239) * q /
        (((((-54.47609879822406 * r + 161.5858368580409) * r - 155.6989798598866) * r + 66.80131188771972) * r - 13.28068155288572) * r + 1);
  }
  for (int j = 0; j < 2; ++j) x += (.5 * std::erfc(x / SQ2) - p) / (std::exp(-.5 * x * x) / std::sqrt(2 * PI));
  return u < .5 ? -x : x;
}


struct Config {
  int logu = 1;              
  double K = 12060.0;        
  double win = 10000.0;      
  double donut_fixed = 500.0;
  double donut_oracle = NaN; 
  double v2 = NaN;           
  int M = 40, R = 499, minpts = 6, run_below = 3, min_run = 20;
  double pthr = 0.10, snr_ceiling = 3.0, se_cap_rel = 0.125, step_u = 0.05;  
  int Kscan = 2000;
  uint64_t seed = 1;         
  double sigma_fixed = NaN;  
  int do_joint_check = 1, do_cpx = 1;
  int reserved_band = 2;         
  int reserved_visibility = 1;           
                             
  double visibility_scale_min = 1.0, u_lo = 0.20, u_hi = 4.00; int minpts_b = 11;
};
enum Status { OK = 0, REF_DENOMINATOR = 1, REF_POOR_FIT = 2, REF_RANGE = 3, REF_INVALID_FIT = 4, REF_VISIBILITY = 5, REF_NOBAND = 6 };
static const char* STATUS_NAME[7] = {"ok", "nonpositive_denominator", "poor_fit", "estimate_outside_range", "fit_not_converged", "no_frequency_band", "no_coherent_window"};
struct Est {
  double est = NaN, p = NaN, se = NaN, den = NaN, pw_mean = NaN, pw_med = NaN, imag = NaN;
  int a = -1, e = -1, status = REF_NOBAND, rank = -1, mask = 0;  
};
struct Result {
  vector<Est> windows;
  size_t n = 0, n_pos = 0;
  double W = NaN, W2 = NaN, ymin = NaN, ymax = NaN, ysd = NaN;
  double sigma_hat = NaN, sig_filt = NaN, sig_iter = NaN; int sigma_fallback = 0;
  double tc = NaN, D = NaN; int nvis = 0, frequency_band_visible = 0, nband = 0, scan_refine = 0;
  double q1 = NaN, q2 = NaN; int g1op = 0, g2op = 0, ng = 0;
  Est cos_raw, cos_chk, cpx_raw, cpx_chk, joint_checked, rat_cos, rat_cpx;
  double pjoint_on_cos = NaN, pass_cos = NaN, pass_cpx = NaN, pass_joint = NaN;
  double pfull_cos = NaN, pfull_cpx = NaN, pfull_joint = NaN;  
  double naive = NaN, donut_fixed = NaN, donut_feas = NaN, donut_oracle = NaN, h_feas = NaN;
  int chol_jitter = 0;
};


static vector<double> least_squares(vector<vector<double>> A, vector<double> y) {  
  int n = A.size(), m = A[0].size();
  for (int k = 0; k < m; k++) {
    double norm = 0; for (int i = k; i < n; i++) norm = std::hypot(norm, A[i][k]); if (norm == 0) return {};
    double alpha = -std::copysign(norm, A[k][k]); vector<double> v(n - k); for (int i = k; i < n; i++) v[i - k] = A[i][k]; v[0] -= alpha;
    double vv = 0; for (double z : v) vv += z * z;
    for (int j = k; j < m; j++) { double q = 0; for (int i = k; i < n; i++) q += v[i - k] * A[i][j]; q *= 2 / vv; for (int i = k; i < n; i++) A[i][j] -= q * v[i - k]; }
    double q = 0; for (int i = k; i < n; i++) q += v[i - k] * y[i]; q *= 2 / vv; for (int i = k; i < n; i++) y[i] -= q * v[i - k];
  }
  vector<double> x(m); for (int i = m - 1; i >= 0; i--) { double z = y[i]; for (int j = i + 1; j < m; j++) z -= A[i][j] * x[j]; x[i] = z / A[i][i]; }
  return x;
}
static bool solve3(const double A[3][3], const double b[3], double x[3]) {  
  double d = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
  double sc = std::fabs(A[0][0] * A[1][1] * A[2][2]);
  if (!(std::fabs(d) > 1e-12 * sc) || !std::isfinite(d)) return false;
  for (int c = 0; c < 3; ++c) {
    double B[3][3]; for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) B[i][j] = (j == c) ? b[i] : A[i][j];
    x[c] = (B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1]) - B[0][1] * (B[1][0] * B[2][2] - B[1][2] * B[2][0]) + B[0][2] * (B[1][0] * B[2][1] - B[1][1] * B[2][0])) / d;
  }
  return std::isfinite(x[0]) && std::isfinite(x[1]) && std::isfinite(x[2]);
}
static bool inv3(const double A[3][3], double I[3][3]) {
  double e[3]; for (int c = 0; c < 3; ++c) { double b[3] = {0, 0, 0}; b[c] = 1; double x[3]; if (!solve3(A, b, x)) return false; for (int r = 0; r < 3; ++r) I[r][c] = x[r]; }
  (void)e; return true;
}
static double median_of(vector<double> x) {
  if (x.empty()) return NaN; std::sort(x.begin(), x.end()); size_t n = x.size();
  return n % 2 ? x[n / 2] : 0.5 * (x[n / 2 - 1] + x[n / 2]);
}





static int rat_rows(bool cplx, const double* cc, const double* ss, const double* c1, const double* s1, const double* cy, const double* sy, int k, double X[2][3], double y[2]) {
  if (!cplx) { X[0][0] = c1[k]; X[0][1] = cc[k] * c1[k]; X[0][2] = -cc[k] * cy[k]; y[0] = cy[k]; return 1; }
  const double zr = cc[k], zi = ss[k];
  X[0][0] = c1[k]; X[0][1] = zr * c1[k] - zi * s1[k]; X[0][2] = -(zr * cy[k] - zi * sy[k]); y[0] = cy[k];
  X[1][0] = s1[k]; X[1][1] = zi * c1[k] + zr * s1[k]; X[1][2] = -(zi * cy[k] + zr * sy[k]); y[1] = sy[k];
  return 2;
}
static bool rat_fit(bool cplx, const double* cc, const double* ss, const double* c1, const double* s1, const double* cy, const double* sy, int a, int e, double th[3]) {
  double A[3][3] = {{0}}, b[3] = {0};
  for (int k = a; k <= e; ++k) { double X[2][3], y[2]; const int nr = rat_rows(cplx, cc, ss, c1, s1, cy, sy, k, X, y);
    for (int q = 0; q < nr; ++q) for (int i = 0; i < 3; ++i) { b[i] += X[q][i] * y[q]; for (int j = 0; j < 3; ++j) A[i][j] += X[q][i] * X[q][j]; } }
  return solve3(A, b, th);
}

static inline double cosU(double sxy, double sxx, double sdY, double sd1) { const double den = sxx - sd1; return den > 0 ? (sxy - sdY) / den : NaN; }
static inline double cpxU(double sam, double saa, int m, double y2, double w2) { const double den = saa - m * w2; return den > 0 ? (sam - m * y2) / den : NaN; }



enum Prod { P11 = 0, P1Y, PYY, P1g1, P1g2, PYg1, PYg2, Pg1g2, NPROD };
static inline int prodidx(int a, int b) {  
  if (a > b) std::swap(a, b);
  if (a == 0) return b == 0 ? P11 : b == 1 ? P1Y : b == 2 ? P1g1 : P1g2;
  if (a == 1) return b == 1 ? PYY : b == 2 ? PYg1 : PYg2;
  if (a == 2) return b == 2 ? P1g1 : Pg1g2;
  return P1g2;  
}
struct Sums {
  int M = 0; double D = 0, W = 0, W2 = 0;
  vector<double> F[4], Fs[2];          
  vector<double> Q[NPROD], Qs[3];      
};

static void main_sums(const vector<double>& V, const vector<double>& om, const vector<double>& Y, const vector<char>& g1, const vector<char>& g2,
                      int ng, double D, int M, Sums& S) {
  const int K2 = 2 * M; const size_t n = V.size();
  S.M = M; S.D = D;
  vector<double> F[4], Fs[2];
  for (int g = 0; g < 4; ++g) F[g].assign(M + 1, 0.0);
  for (int g = 0; g < 2; ++g) Fs[g].assign(M + 1, 0.0);
  for (int p = 0; p < NPROD; ++p) S.Q[p].assign(K2 + 1, 0.0);
  for (int p = 0; p < 3; ++p) S.Qs[p].assign(K2 + 1, 0.0);
  constexpr int B = 256;
  alignas(64) double c[B], s[B], cd[B], sd[B], w1[4][B], w2[NPROD][B];
  for (size_t i0 = 0; i0 < n; i0 += B) {
    const int L = (int)std::min<size_t>(B, n - i0);
    for (int ii = 0; ii < L; ++ii) {
      const size_t i = i0 + ii; const double o = om[i], y = Y[i], a1 = (ng > 0 && g1[i]) ? 1.0 : 0.0, a2 = (ng > 1 && g2[i]) ? 1.0 : 0.0, o2 = o * o;
      const double x = D * V[i]; cd[ii] = std::cos(x); sd[ii] = std::sin(x); c[ii] = 1.0; s[ii] = 0.0;
      w1[0][ii] = o; w1[1][ii] = o * y; w1[2][ii] = o * a1; w1[3][ii] = o * a2;
      w2[P11][ii] = o2; w2[P1Y][ii] = o2 * y; w2[PYY][ii] = o2 * y * y; w2[P1g1][ii] = o2 * a1; w2[P1g2][ii] = o2 * a2;
      w2[PYg1][ii] = o2 * y * a1; w2[PYg2][ii] = o2 * y * a2; w2[Pg1g2][ii] = o2 * a1 * a2;
    }
    const int nfo = 2 + ng, nso = (ng == 0) ? 3 : (ng == 1 ? 5 : NPROD);
    static const int so_list[NPROD] = {P11, P1Y, PYY, P1g1, PYg1, P1g2, PYg2, Pg1g2};
    for (int k = 0; k <= K2; ++k) {
      if (k <= M) {
        for (int g = 0; g < nfo; ++g) {
          const double* w = w1[g]; double acc = 0;
#pragma omp simd reduction(+:acc)
          for (int ii = 0; ii < L; ++ii) acc += w[ii] * c[ii];
          F[g][k] += acc;
        }
        for (int g = 0; g < 2; ++g) {
          const double* w = w1[g]; double acc = 0;
#pragma omp simd reduction(+:acc)
          for (int ii = 0; ii < L; ++ii) acc += w[ii] * s[ii];
          Fs[g][k] += acc;
        }
      }
      for (int q = 0; q < nso; ++q) {
        const int p = so_list[q]; const double* w = w2[p]; double acc = 0;
#pragma omp simd reduction(+:acc)
        for (int ii = 0; ii < L; ++ii) acc += w[ii] * c[ii];
        S.Q[p][k] += acc;
      }
      for (int p = 0; p < 3; ++p) {
        const double* w = w2[p]; double acc = 0;
#pragma omp simd reduction(+:acc)
        for (int ii = 0; ii < L; ++ii) acc += w[ii] * s[ii];
        S.Qs[p][k] += acc;
      }
#pragma omp simd
      for (int ii = 0; ii < L; ++ii) { const double cn = c[ii] * cd[ii] - s[ii] * sd[ii], sn = s[ii] * cd[ii] + c[ii] * sd[ii]; c[ii] = cn; s[ii] = sn; }
    }
  }
  S.W = F[0][0]; S.W2 = S.Q[P11][0];
  for (int g = 0; g < 4; ++g) { S.F[g].assign(M + 1, 0.0); for (int k = 0; k <= M; ++k) S.F[g][k] = F[g][k] / S.W; }
  for (int g = 0; g < 2; ++g) { S.Fs[g].assign(M + 1, 0.0); for (int k = 0; k <= M; ++k) S.Fs[g][k] = Fs[g][k] / S.W; }
}




static int ceiling_scan(const vector<double>& V, const vector<double>& om, double D0, int Kscan, double thr, int run_below, int& scanned) {
  const size_t n = V.size(); constexpr int CH = 64;
  vector<double> c(n, 1.0), s(n, 0.0), cd(n), sd(n);
  double W = 0, W2 = 0;
  for (size_t i = 0; i < n; ++i) { cd[i] = std::cos(D0 * V[i]); sd[i] = std::sin(D0 * V[i]); W += om[i]; W2 += om[i] * om[i]; }
  int run = 0, last_vis = 0; scanned = 0;
  for (int k0 = 1; k0 <= Kscan; k0 += CH) {
    const int kk = std::min(CH, Kscan - k0 + 1);
    double A1[CH] = {0}, A2[CH] = {0}, B1[CH] = {0}, B2[CH] = {0};
    for (size_t i0 = 0; i0 < n; i0 += 256) {
      const int L = (int)std::min<size_t>(256, n - i0);
      double* cc = &c[i0]; double* ss = &s[i0]; const double* cdd = &cd[i0]; const double* sdd = &sd[i0]; const double* o = &om[i0];
      for (int j = 0; j < kk; ++j) {
        double a1 = 0, a2 = 0, b1 = 0, b2 = 0;
#pragma omp simd reduction(+:a1,a2,b1,b2)
        for (int ii = 0; ii < L; ++ii) {
          const double cn = cc[ii] * cdd[ii] - ss[ii] * sdd[ii], sn = ss[ii] * cdd[ii] + cc[ii] * sdd[ii];
          cc[ii] = cn; ss[ii] = sn;
          const double oc = o[ii] * cn, os = o[ii] * sn; a1 += oc; a2 += o[ii] * oc; b1 += os; b2 += o[ii] * os;
        }
        A1[j] += a1; A2[j] += a2; B1[j] += b1; B2[j] += b2;
      }
    }
    for (int j = 0; j < kk; ++j) {
      const int k = k0 + j; ++scanned;
      const double C = A1[j] / W, Sn = B1[j] / W, a2 = C * C + Sn * Sn;
      const double var = (W2 * (1 + a2) - 2 * (C * A2[j] + Sn * B2[j])) / (W * W);
      const double snr = std::sqrt(a2) / std::sqrt(std::max(var, 1e-300));
      if (snr >= thr) { run = 0; last_vis = k; }
      else if (++run >= run_below) return last_vis;
    }
  }
  return last_vis;
}


struct Comp { int g, trig; };  
struct Cov {
  const Sums* S; vector<Comp> comps; int M, N; vector<double> A;  
  double cosv(int p, int j) const { return S->Q[p][std::abs(j)]; }
  double sinv(int p, int j) const { return j >= 0 ? S->Qs[p][j] : -S->Qs[p][-j]; }
  double X(const Comp& c, int k) const { return c.trig == 0 ? S->F[c.g][k] : S->Fs[c.g][k]; }
  double F(const Comp& c, int k) const { int p = prodidx(c.g, 0); return c.trig == 0 ? cosv(p, k) : sinv(p, k); }
  double cov(const Comp& p, int k, const Comp& q, int l) const {
    const int pr = prodidx(p.g, q.g); double E;
    if (p.trig == 0 && q.trig == 0) E = 0.5 * (cosv(pr, k - l) + cosv(pr, k + l));
    else if (p.trig == 1 && q.trig == 1) E = 0.5 * (cosv(pr, k - l) - cosv(pr, k + l));
    else if (p.trig == 0) E = 0.5 * (sinv(pr, l + k) + sinv(pr, l - k));
    else E = 0.5 * (sinv(pr, k + l) + sinv(pr, k - l));
    const double Xp = X(p, k), Xq = X(q, l);
    return (E - Xq * F(p, k) - Xp * F(q, l) + Xp * Xq * S->W2) / (S->W * S->W);
  }
  void build(const Sums& s, const vector<Comp>& cs) {
    S = &s; comps = cs; M = s.M; N = (int)cs.size() * M; A.assign((size_t)N * N, 0.0);
    for (int a = 0; a < (int)cs.size(); ++a) for (int k = 1; k <= M; ++k)
      for (int b = 0; b < (int)cs.size(); ++b) for (int l = 1; l <= M; ++l)
        A[(size_t)(a * M + k - 1) * N + (b * M + l - 1)] = cov(cs[a], k, cs[b], l);
  }
  double at(int a, int k, int b, int l) const { return A[(size_t)(a * M + k - 1) * N + (b * M + l - 1)]; }  
};

static int cholesky(const vector<double>& A, int N, vector<double>& L) {
  double tr = 0; for (int i = 0; i < N; ++i) tr += A[(size_t)i * N + i];
  double jit = 1e-10 * tr / N;
  for (int attempt = 0; attempt < 4; ++attempt, jit *= 100) {
    L.assign((size_t)N * N, 0.0); bool ok = true;
    for (int i = 0; i < N && ok; ++i) for (int j = 0; j <= i; ++j) {
      double s = A[(size_t)i * N + j];
      for (int k = 0; k < j; ++k) s -= L[(size_t)i * N + k] * L[(size_t)j * N + k];
      if (i == j) { s += jit; if (!(s > 0)) { ok = false; break; } L[(size_t)i * N + i] = std::sqrt(s); }
      else L[(size_t)i * N + j] = s / L[(size_t)j * N + j];
    }
    if (ok) return attempt;
  }
  return -1;
}


struct Band { int a, e, m; bool valid = false; double mu = NaN, est = NaN, den = NaN, p = NaN, pj = NaN; };
static void enumerate_bands(int M, int minpts, vector<Band>& B, int amin = 0) {
  B.clear();
  for (int a = amin; a + minpts <= M; ++a) for (int e = a + minpts - 1; e < M; ++e) { Band b; b.a = a; b.e = e; b.m = e - a + 1; B.push_back(b); }
}



template <class P, class DenominatorCheck, class RangeCheck>
static void select_ranked_windows(const vector<Band>& B, P p_of, DenominatorCheck positive_denominator, RangeCheck within_range, double pthr, bool frequency_band_visible, Est& raw, Est& chk) {
  vector<int> hi, lo;
  for (int i = 0; i < (int)B.size(); ++i) if (B[i].valid) (p_of(B[i]) >= pthr ? hi : lo).push_back(i);
  std::stable_sort(hi.begin(), hi.end(), [&](int x, int y) { if (B[x].m != B[y].m) return B[x].m > B[y].m; return B[x].a < B[y].a; });
  std::stable_sort(lo.begin(), lo.end(), [&](int x, int y) { double px = p_of(B[x]), py = p_of(B[y]); if (px != py) return px > py; if (B[x].m != B[y].m) return B[x].m > B[y].m; return B[x].a < B[y].a; });
  vector<int> rk(hi); rk.insert(rk.end(), lo.begin(), lo.end());
  auto fill = [&](Est& o, int i, int rank) { const Band& b = B[i]; o.a = b.a; o.e = b.e; o.p = p_of(b); o.den = b.den; o.est = b.est; o.rank = rank; };
  raw = Est(); chk = Est();
  if (!rk.empty()) { fill(raw, rk[0], 1); raw.status = OK; }
  int mask = frequency_band_visible ? 1 : 0;
  bool s1 = false, s2 = false, s3 = false; int ans = -1, ansrank = -1;
  for (int i : rk) if (positive_denominator(B[i])) { s1 = true; break; }
  for (int i : hi) if (positive_denominator(B[i])) { s2 = true; break; }
  for (int r = 0; r < (int)hi.size(); ++r) {
    const int i = hi[r];
    if (!positive_denominator(B[i]) || !within_range(B[i])) continue;
    s3 = true;
    ans = i; ansrank = r + 1; break;
  }
  if (s1) mask |= 2; if (s2) mask |= 4; if (s3) mask |= 8; if (ans >= 0) mask |= 16;
  chk.mask = mask; raw.mask = mask;
  if (!frequency_band_visible) { chk.status = REF_VISIBILITY; raw.status = rk.empty() ? REF_NOBAND : raw.status; return; }
  if (rk.empty()) { chk.status = REF_NOBAND; return; }
  if (ans >= 0) { fill(chk, ans, ansrank); chk.status = OK; return; }
  chk.status = !s1 ? REF_DENOMINATOR : !s2 ? REF_POOR_FIT : REF_RANGE;
}


static double donut_h(const vector<double>& V, const vector<double>& Y, const vector<double>& cw, double sc, double h, double win) {
  
  
  
  if (!(h > 0) || !std::isfinite(h) || !std::isfinite(win) || !(win > 0) || !std::isfinite(sc)) return NaN;
  const double Hraw = win + (h - 500.0); const int nb = 2 * (int)std::ceil(Hraw / 250.0 - 1e-9); const double H = nb * 125.0;
  vector<double> hc(nb, 0.0), hy(nb, 0.0);
  for (size_t i = 0; i < V.size(); ++i) {
    if (!(cw[i] > 0)) continue;
    const int hb = int(std::floor((V[i] * sc + H) / 250));
    if (hb >= 0 && hb < nb) { hc[hb] += cw[i]; hy[hb] += cw[i] * Y[i]; }
  }
  vector<vector<double>> A; vector<double> yy;
  for (int i = 0; i < nb; i++) { double x = (i + .5) * 250 - H; if (std::abs(x) <= h) continue; vector<double> row(8, 1); for (int j = 1; j < 8; j++) row[j] = row[j - 1] * x / H; A.push_back(row); yy.push_back(hc[i]); }
  if ((int)A.size() < 8) return NaN;
  auto coef = least_squares(A, yy); if (coef.empty()) return NaN;
  const double side = h + 2500.0;
  vector<double> lines[2];
  for (int sd = 0; sd < 2; sd++) {
    A.clear(); yy.clear();
    for (int i = 0; i < nb; i++) { double x = (i + .5) * 250 - H; if (std::abs(x) <= h || std::abs(x) > side || (x > 0) != (sd == 1) || hc[i] == 0) continue; double w = std::sqrt(hc[i]); A.push_back({w, w * x / H}); yy.push_back(hy[i] / hc[i] * w); }
    if (A.size() < 2) return NaN;
    lines[sd] = least_squares(A, yy); if (lines[sd].empty()) return NaN;
  }
  double num = 0, den = 0;
  for (int i = 0; i < nb; i++) { double x = (i + .5) * 250 - H; if (std::abs(x) > h) continue; double z = 0; for (int j = 7; j >= 0; j--) z = z * x / H + coef[j]; auto& l = lines[x > 0]; num += hy[i] - z * (l[0] + l[1] * x / H); den += hc[i] - z; }
  return den > 0 ? num / den : NaN;
}




static void cell_weights(size_t n, const double* w, const double* tau, int weighted, vector<double>& om, vector<double>& cw) {
  om.resize(n); cw.resize(n);
  for (size_t i = 0; i < n; ++i) { const double wi = (weighted && w) ? w[i] : 1.0; om[i] = wi * tau[i]; cw[i] = tau[i] > 0 ? wi : 0.0; }
}


static Result analyze(const vector<double>& Vall, const vector<double>& Yall, const vector<double>& omall, const vector<double>& cwall, const Config& cf) {
  Result R; R.n = Vall.size();
  const double sc = cf.logu ? cf.K : 1.0;  
  
  {
    double sw = 0, sy = 0; for (size_t i = 0; i < Vall.size(); ++i) if (cwall[i] > 0 && std::fabs(Vall[i] * sc) <= cf.win) { sw += cwall[i]; sy += cwall[i] * Yall[i]; }
    R.naive = sw > 0 ? sy / sw : NaN;
    R.donut_fixed = donut_h(Vall, Yall, cwall, sc, cf.donut_fixed, cf.win);
    R.donut_oracle = std::isfinite(cf.donut_oracle) ? donut_h(Vall, Yall, cwall, sc, cf.donut_oracle, cf.win) : NaN;
  }
  
  vector<double> V, Y, om; V.reserve(Vall.size()); Y.reserve(Vall.size()); om.reserve(Vall.size());
  for (size_t i = 0; i < Vall.size(); ++i) if (omall[i] > 0) { V.push_back(Vall[i]); Y.push_back(Yall[i]); om.push_back(omall[i]); }
  R.n_pos = V.size();
  if (V.size() < 100) return R;
  {
    double sw = 0, sy = 0, syy = 0; R.ymin = INFINITY; R.ymax = -INFINITY;
    for (size_t i = 0; i < V.size(); ++i) { sw += om[i]; sy += om[i] * Y[i]; R.ymin = std::min(R.ymin, Y[i]); R.ymax = std::max(R.ymax, Y[i]); }
    const double mY = sy / sw; for (size_t i = 0; i < V.size(); ++i) syy += om[i] * (Y[i] - mY) * (Y[i] - mY);
    R.ysd = std::sqrt(syy / sw);
  }
  
  if (std::isfinite(cf.sigma_fixed)) R.sigma_hat = cf.sigma_fixed;
  else {
    vector<double> dZ(V.size()); const double f = cf.logu ? 1.0 : 1.0 / cf.K;
    for (size_t i = 0; i < V.size(); ++i) dZ[i] = V[i] * f;
    p2pilot::Pilot P = p2pilot::ensemble_sigma(dZ);
    R.sig_filt = P.sig_filt; R.sig_iter = P.sig_iter;
    R.sigma_hat = std::isfinite(P.sigma_p) ? P.sigma_p / f : NaN;
  }
  double sig = R.sigma_hat;
  if (!(std::isfinite(sig) && sig > 0)) { sig = 1000.0 / sc; R.sigma_fallback = 1; }  
  R.h_feas = std::isfinite(R.sigma_hat) ? 2.0 * R.sigma_hat * sc : NaN;
  R.donut_feas = donut_h(Vall, Yall, cwall, sc, R.h_feas, cf.win);
  
  double D0 = cf.step_u / sig; int scanned = 0, cidx = 0;
  for (int ref = 0; ref <= 4; ++ref) {
    cidx = ceiling_scan(V, om, D0, cf.Kscan, cf.snr_ceiling, cf.run_below, scanned);
    R.scan_refine = ref;
    if (cidx >= cf.min_run || ref == 4) break;
    D0 *= 0.5;
  }
  R.nvis = cidx; R.tc = cidx * D0;
  R.frequency_band_visible = cf.reserved_visibility == 0 ? (cidx >= 2) : cf.reserved_visibility == 1 ? 1 : (sig * R.tc >= cf.visibility_scale_min);
  
  const bool part1grid = cf.reserved_band == 2;
  const int M = part1grid ? (int)std::lround(cf.u_hi / cf.step_u) : cf.M;
  const double D = part1grid ? cf.step_u / sig : R.tc / M;
  const int amin = part1grid ? (int)std::lround(cf.u_lo / cf.step_u) - 1 : 0, minpts = part1grid ? cf.minpts_b : cf.minpts;
  const bool gridok = part1grid || cidx >= minpts;  
  if (!R.frequency_band_visible || !gridok) {
    const int st = !R.frequency_band_visible ? REF_VISIBILITY : REF_NOBAND;
    R.cos_chk.status = R.cpx_chk.status = R.joint_checked.status = R.rat_cos.status = R.rat_cpx.status = st;
    return R;
  }
  
  vector<char> g1(V.size(), 0), g2(V.size(), 0); int ng = 0;
  if (cf.do_joint_check) {
    vector<std::pair<double, double>> yw; double tot = 0;
    for (size_t i = 0; i < V.size(); ++i) if (std::fabs(V[i]) <= 3 * sig) { yw.push_back({Y[i], om[i]}); tot += om[i]; }
    if (yw.size() < 50) { for (size_t i = 0; i < V.size(); ++i) yw.push_back({Y[i], om[i]}); tot = 0; for (auto& q : yw) tot += q.second; }
    std::sort(yw.begin(), yw.end());
    auto wq = [&](double p) { double cum = 0; for (auto& q : yw) { cum += q.second; if (cum >= p * tot) return q.first; } return yw.back().first; };
    const double qs[2] = {wq(1.0 / 3), wq(2.0 / 3)};
    double wt = 0; for (size_t i = 0; i < V.size(); ++i) wt += om[i];
    vector<char>* slots[2] = {&g1, &g2}; double prev_share = NaN;
    for (int j = 0; j < 2; ++j) {
      int chosen = -1; double share = NaN;
      for (int op = 0; op < 2 && chosen < 0; ++op) {  
        double s1 = 0; for (size_t i = 0; i < V.size(); ++i) if (op == 0 ? Y[i] <= qs[j] : Y[i] < qs[j]) s1 += om[i];
        const double sh = s1 / wt; if (sh > 1e-3 && sh < 1 - 1e-3) { chosen = op; share = sh; }
      }
      if (chosen < 0) continue;                                            
      if (ng == 1 && std::fabs(share - prev_share) < 1e-12) continue;      
      vector<char>& g = *slots[ng];
      for (size_t i = 0; i < V.size(); ++i) g[i] = (chosen == 0 ? Y[i] <= qs[j] : Y[i] < qs[j]) ? 1 : 0;
      if (ng == 0) { R.q1 = qs[j]; R.g1op = chosen; } else { R.q2 = qs[j]; R.g2op = chosen; }
      prev_share = share; ++ng;
    }
  }
  R.ng = ng;
  
  R.D = D;
  Sums S; main_sums(V, om, Y, g1, g2, ng, D, M, S); R.W = S.W; R.W2 = S.W2;
  
  
  enum { iC1 = 0, iCY = 1, iS1 = 2, iSY = 3, iG1 = 2, iG2 = 3 };
  const int Rb = cf.R;
  struct Draws { Cov CV; vector<double> Z; int jit = -1; int N = 0; };
  auto make = [&](const vector<Comp>& comps, uint64_t salt, Draws& d) {
    d.CV.build(S, comps); d.N = d.CV.N; vector<double> L; d.jit = cholesky(d.CV.A, d.N, L);
    d.Z.assign((size_t)Rb * d.N, 0.0);
    if (d.jit < 0) return;
    vector<double> xi(d.N);
    for (int r = 0; r < Rb; ++r) {
      CStream st(sm64(cf.seed ^ salt), (uint64_t)r);
      for (int i = 0; i < d.N; ++i) xi[i] = Qn(st.u());
      double* z = &d.Z[(size_t)r * d.N];
      for (int i = 0; i < d.N; ++i) { double acc = 0; const double* Li = &L[(size_t)i * d.N]; for (int k = 0; k <= i; ++k) acc += Li[k] * xi[k]; z[i] = acc; }
    }
  };
  Draws DA, DB;
  { vector<Comp> ca = {{0, 0}, {1, 0}}; for (int j = 0; j < ng; ++j) ca.push_back({2 + j, 0}); make(ca, 0xB00757A7ULL, DA); }
  if (cf.do_cpx || cf.reserved_band == 1 || std::isfinite(cf.v2)) make({{0, 0}, {1, 0}, {0, 1}, {1, 1}}, 0xB00757B8ULL, DB);
  R.chol_jitter = DA.jit;
  auto zA = [&](int r, int comp, int k) { return DA.Z[(size_t)r * DA.N + comp * M + k]; };  
  auto zB = [&](int r, int comp, int k) { return DB.Z[(size_t)r * DB.N + comp * M + k]; };
  
  vector<double> cosine_mean(M), CY(M), S1(M), SY(M), CG[2];
  for (int k = 0; k < M; ++k) { cosine_mean[k] = S.F[0][k + 1]; CY[k] = S.F[1][k + 1]; S1[k] = S.Fs[0][k + 1]; SY[k] = S.Fs[1][k + 1]; }
  for (int j = 0; j < ng; ++j) { CG[j].resize(M); for (int k = 0; k < M; ++k) CG[j][k] = S.F[2 + j][k + 1]; }
  const double WW = S.W * S.W;
  vector<double> d1(M), dY(M);
  for (int k = 0; k < M; ++k) { d1[k] = (S.Q[P11][0] + S.Q[P11][2 * (k + 1)]) / (2 * WW); dY[k] = (S.Q[P1Y][0] + S.Q[P1Y][2 * (k + 1)]) / (2 * WW); }
  auto VA = [&](int a, int k, int b, int l) { return DA.CV.at(a, k + 1, b, l + 1); };
  auto VB = [&](int a, int k, int b, int l) { return DB.CV.at(a, k + 1, b, l + 1); };
  
  auto prefix = [&](auto f) { vector<double> p(M + 1, 0.0); for (int k = 0; k < M; ++k) p[k + 1] = p[k] + f(k); return p; };
  const vector<double> pXY = prefix([&](int k) { return cosine_mean[k] * CY[k]; }), pXX = prefix([&](int k) { return cosine_mean[k] * cosine_mean[k]; });
  const vector<double> pD1 = prefix([&](int k) { return d1[k]; }), pDY = prefix([&](int k) { return dY[k]; });
  const vector<double> pAM = prefix([&](int k) { return cosine_mean[k] * CY[k] + S1[k] * SY[k]; }), pAA = prefix([&](int k) { return cosine_mean[k] * cosine_mean[k] + S1[k] * S1[k]; });
  vector<Band> B; enumerate_bands(M, minpts, B, amin); R.nband = (int)B.size();
  const bool drawA = DA.jit >= 0, drawB = DB.jit >= 0;
  
  vector<Band> Bc = B;
  {
    vector<double> Pa((size_t)Rb * (M + 1), 0.0), Pb(Pa), Pc(Pa);
    for (int r = 0; r < Rb; ++r) for (int k = 0; k < M; ++k) {
      const double c1s = cosine_mean[k] + zA(r, iC1, k); const size_t o = (size_t)r * (M + 1) + k;
      Pa[o + 1] = Pa[o] + c1s * cosine_mean[k]; Pb[o + 1] = Pb[o] + c1s * zA(r, iCY, k); Pc[o + 1] = Pc[o] + c1s * c1s;
    }
    vector<double> sdv(M);
    int npass = 0, nval = 0;
    for (auto& b : Bc) {
      const double sxx = pXX[b.e + 1] - pXX[b.a], sxy = pXY[b.e + 1] - pXY[b.a];
      if (!(sxx > 0)) continue;
      const double mu = sxy / sxx; if (!std::isfinite(mu)) continue;
      b.valid = true; b.mu = mu; ++nval;
      b.den = sxx - (pD1[b.e + 1] - pD1[b.a]);
      b.est = cosU(sxy, sxx, pDY[b.e + 1] - pDY[b.a], pD1[b.e + 1] - pD1[b.a]);
      double T = 0;
      for (int k = b.a; k <= b.e; ++k) {
        const double v = VA(iCY, k, iCY, k) - 2 * mu * VA(iCY, k, iC1, k) + mu * mu * VA(iC1, k, iC1, k);
        sdv[k] = std::sqrt(std::max(v, 1e-300)); T = std::max(T, std::fabs(CY[k] - mu * cosine_mean[k]) / sdv[k]);
      }
      int ex = 0;
      if (drawA) for (int r = 0; r < Rb; ++r) {
        const size_t o = (size_t)r * (M + 1);
        const double ms = (mu * (Pa[o + b.e + 1] - Pa[o + b.a]) + (Pb[o + b.e + 1] - Pb[o + b.a])) / (Pc[o + b.e + 1] - Pc[o + b.a]);
        if (!std::isfinite(ms)) continue;
        for (int k = b.a; k <= b.e; ++k) {
          const double c1s = cosine_mean[k] + zA(r, iC1, k), cys = mu * cosine_mean[k] + zA(r, iCY, k);
          if (std::fabs(cys - ms * c1s) / sdv[k] >= T) { ++ex; break; }
        }
      }
      b.p = drawA ? (1.0 + ex) / (Rb + 1.0) : NaN;
      if (b.p >= cf.pthr) ++npass;
    }
    R.pass_cos = nval ? (double)npass / nval : NaN;
  }
  
  if (cf.do_joint_check) {
    const int G = 1 + ng; const int gcomp[3] = {iCY, iG1, iG2}; const vector<double>* Cg[3] = {&CY, &CG[0], &CG[1]};
    vector<double> Pa((size_t)Rb * (M + 1), 0.0), Pc(Pa), Pb[3];
    for (int g = 0; g < G; ++g) Pb[g].assign((size_t)Rb * (M + 1), 0.0);
    for (int r = 0; r < Rb; ++r) for (int k = 0; k < M; ++k) {
      const double c1s = cosine_mean[k] + zA(r, iC1, k); const size_t o = (size_t)r * (M + 1) + k;
      Pa[o + 1] = Pa[o] + c1s * cosine_mean[k]; Pc[o + 1] = Pc[o] + c1s * c1s;
      for (int g = 0; g < G; ++g) Pb[g][o + 1] = Pb[g][o] + c1s * zA(r, gcomp[g], k);
    }
    vector<vector<double>> pg(G);
    for (int g = 0; g < G; ++g) { const vector<double>& cg = *Cg[g]; pg[g] = prefix([&](int k) { return cosine_mean[k] * cg[k]; }); }
    vector<double> sdv((size_t)G * M);
    int npass = 0, nval = 0;
    for (auto& b : Bc) {
      if (!b.valid) continue;
      const double sxx = pXX[b.e + 1] - pXX[b.a]; double mug[3]; double T = 0; ++nval;
      for (int g = 0; g < G; ++g) {
        const vector<double>& cg = *Cg[g]; mug[g] = (pg[g][b.e + 1] - pg[g][b.a]) / sxx;
        for (int k = b.a; k <= b.e; ++k) {
          const double v = VA(gcomp[g], k, gcomp[g], k) - 2 * mug[g] * VA(gcomp[g], k, iC1, k) + mug[g] * mug[g] * VA(iC1, k, iC1, k);
          sdv[(size_t)g * M + k] = std::sqrt(std::max(v, 1e-300)); T = std::max(T, std::fabs(cg[k] - mug[g] * cosine_mean[k]) / sdv[(size_t)g * M + k]);
        }
      }
      int ex = 0;
      if (drawA) for (int r = 0; r < Rb; ++r) {
        const size_t o = (size_t)r * (M + 1); const double dA = Pa[o + b.e + 1] - Pa[o + b.a], dC = Pc[o + b.e + 1] - Pc[o + b.a];
        bool hit = false;
        for (int g = 0; g < G && !hit; ++g) {
          const double ms = (mug[g] * dA + (Pb[g][o + b.e + 1] - Pb[g][o + b.a])) / dC;
          if (!std::isfinite(ms)) continue;
          for (int k = b.a; k <= b.e; ++k) {
            const double c1s = cosine_mean[k] + zA(r, iC1, k), cgs = mug[g] * cosine_mean[k] + zA(r, gcomp[g], k);
            if (std::fabs(cgs - ms * c1s) / sdv[(size_t)g * M + k] >= T) { hit = true; break; }
          }
        }
        ex += hit;
      }
      b.pj = drawA ? (1.0 + ex) / (Rb + 1.0) : NaN;
      if (b.pj >= cf.pthr) ++npass;
    }
    R.pass_joint = nval ? (double)npass / nval : NaN;
  }
  
  auto se_cos = [&](const Band& b) {
    const double mu = b.mu, sxx = pXX[b.e + 1] - pXX[b.a]; double v = 0;
    for (int k = b.a; k <= b.e; ++k) for (int l = b.a; l <= b.e; ++l) {
      const double om_ = VA(iCY, k, iCY, l) - mu * VA(iCY, k, iC1, l) - mu * VA(iC1, k, iCY, l) + mu * mu * VA(iC1, k, iC1, l);
      v += cosine_mean[k] * cosine_mean[l] * om_;
    }
    return std::sqrt(std::max(v, 0.0)) / sxx;
  };
  auto positive_denominator = [&](const Band& b) { return b.den > 0 && std::isfinite(b.est); };
  auto within_range = [&](const Band& b) { return b.est >= R.ymin && b.est <= R.ymax; };
  auto pw = [&](Est& o) {
    if (o.a < 0) return; vector<double> r; double im = 0;
    for (int k = o.a; k <= o.e; ++k) { r.push_back(CY[k] / cosine_mean[k]); const double a2 = cosine_mean[k] * cosine_mean[k] + S1[k] * S1[k]; im += (SY[k] * cosine_mean[k] - CY[k] * S1[k]) / a2; }
    double s = 0; for (double x : r) s += x; o.pw_mean = s / r.size(); o.pw_med = median_of(r); o.imag = im / r.size();
  };
  auto fin_cos = [&](Est& o) { if (o.a >= 0) { for (auto& b : Bc) if (b.a == o.a && b.e == o.e) { o.se = se_cos(b); break; } } pw(o); };
  select_ranked_windows(Bc, [](const Band& b) { return b.p; }, positive_denominator, within_range, cf.pthr, R.frequency_band_visible, R.cos_raw, R.cos_chk);
  fin_cos(R.cos_raw); fin_cos(R.cos_chk);
  vector<int> ranked;
  for (int i=0;i<(int)Bc.size();++i) if(Bc[i].valid) ranked.push_back(i);
  std::stable_sort(ranked.begin(),ranked.end(),[&](int i,int j) {
    const auto &a=Bc[i], &b=Bc[j];
    bool ap=a.p>=cf.pthr, bp=b.p>=cf.pthr;
    if(ap!=bp) return ap;
    if(!ap && a.p!=b.p) return a.p>b.p;
    if(a.m!=b.m) return a.m>b.m;
    return a.a<b.a;
  });
  for(int i:ranked) {
    const auto &b=Bc[i]; Est e;
    e.a=b.a;e.e=b.e;e.p=b.p;e.den=b.den;e.est=b.est;e.rank=(int)R.windows.size()+1;
    e.status=positive_denominator(b)&&within_range(b)?OK:REF_DENOMINATOR;
    fin_cos(e);R.windows.push_back(e);
  }
  if(R.cos_chk.status!=OK) for(const auto &e:R.windows) if(e.status==OK) {R.cos_chk=e;break;}

  if (R.cos_chk.a >= 0) for (auto& b : Bc) if (b.a == R.cos_chk.a && b.e == R.cos_chk.e) { R.pjoint_on_cos = b.pj; break; }
  for (auto& b : Bc) if (b.a == amin && b.e == M - 1) { R.pfull_cos = b.p; R.pfull_joint = b.pj; break; }
  if (cf.do_joint_check) { Est dummy; select_ranked_windows(Bc, [](const Band& b) { return b.pj; }, positive_denominator, within_range, cf.pthr, R.frequency_band_visible, dummy, R.joint_checked); fin_cos(R.joint_checked); }
  
  if (cf.do_cpx || cf.reserved_band == 1) {
    vector<Band> Bx = B;
    vector<double> Pa((size_t)Rb * (M + 1), 0.0), Pb(Pa), Pc(Pa);
    for (int r = 0; r < Rb; ++r) for (int k = 0; k < M; ++k) {
      const double c1s = cosine_mean[k] + zB(r, iC1, k), s1s = S1[k] + zB(r, iS1, k); const size_t o = (size_t)r * (M + 1) + k;
      Pa[o + 1] = Pa[o] + c1s * cosine_mean[k] + s1s * S1[k]; Pb[o + 1] = Pb[o] + c1s * zB(r, iCY, k) + s1s * zB(r, iSY, k); Pc[o + 1] = Pc[o] + c1s * c1s + s1s * s1s;
    }
    vector<double> sdr(M), sdi(M);
    const double y2 = S.Q[P1Y][0] / WW, w2 = S.W2 / WW;
    int npass = 0, nval = 0;
    for (auto& b : Bx) {
      const double saa = pAA[b.e + 1] - pAA[b.a], sam = pAM[b.e + 1] - pAM[b.a];
      if (!(saa > 0)) continue;
      const double mu = sam / saa; if (!std::isfinite(mu)) continue;
      b.valid = true; b.mu = mu; ++nval;
      b.den = saa - b.m * w2; b.est = cpxU(sam, saa, b.m, y2, w2);
      double T = 0;
      for (int k = b.a; k <= b.e; ++k) {
        const double vr = VB(iCY, k, iCY, k) - 2 * mu * VB(iCY, k, iC1, k) + mu * mu * VB(iC1, k, iC1, k);
        const double vi = VB(iSY, k, iSY, k) - 2 * mu * VB(iSY, k, iS1, k) + mu * mu * VB(iS1, k, iS1, k);
        sdr[k] = std::sqrt(std::max(vr, 1e-300)); sdi[k] = std::sqrt(std::max(vi, 1e-300));
        T = std::max(T, std::max(std::fabs(CY[k] - mu * cosine_mean[k]) / sdr[k], std::fabs(SY[k] - mu * S1[k]) / sdi[k]));
      }
      int ex = 0;
      if (drawB) for (int r = 0; r < Rb; ++r) {
        const size_t o = (size_t)r * (M + 1);
        const double ms = (mu * (Pa[o + b.e + 1] - Pa[o + b.a]) + (Pb[o + b.e + 1] - Pb[o + b.a])) / (Pc[o + b.e + 1] - Pc[o + b.a]);
        if (!std::isfinite(ms)) continue;
        for (int k = b.a; k <= b.e; ++k) {
          const double c1s = cosine_mean[k] + zB(r, iC1, k), s1s = S1[k] + zB(r, iS1, k), cys = mu * cosine_mean[k] + zB(r, iCY, k), sys = mu * S1[k] + zB(r, iSY, k);
          if (std::fabs(cys - ms * c1s) / sdr[k] >= T || std::fabs(sys - ms * s1s) / sdi[k] >= T) { ++ex; break; }
        }
      }
      b.p = drawB ? (1.0 + ex) / (Rb + 1.0) : NaN;
      if (b.p >= cf.pthr) ++npass;
    }
    R.pass_cpx = nval ? (double)npass / nval : NaN;
    for (auto& b : Bx) if (b.a == amin && b.e == M - 1) { R.pfull_cpx = b.p; break; }
    auto se_cpx = [&](const Band& b) {
      const double mu = b.mu, saa = pAA[b.e + 1] - pAA[b.a]; double v = 0;
      
      const int cy[2] = {iCY, iSY}, c1[2] = {iC1, iS1}; const vector<double>* av[2] = {&cosine_mean, &S1};
      for (int u = 0; u < 2; ++u) for (int w = 0; w < 2; ++w)
        for (int k = b.a; k <= b.e; ++k) for (int l = b.a; l <= b.e; ++l) {
          const double om_ = VB(cy[u], k, cy[w], l) - mu * VB(cy[u], k, c1[w], l) - mu * VB(c1[u], k, cy[w], l) + mu * mu * VB(c1[u], k, c1[w], l);
          v += (*av[u])[k] * (*av[w])[l] * om_;
        }
      return std::sqrt(std::max(v, 0.0)) / saa;
    };
    select_ranked_windows(Bx, [](const Band& b) { return b.p; }, positive_denominator, within_range, cf.pthr, R.frequency_band_visible, R.cpx_raw, R.cpx_chk);
    for (Est* o : {&R.cpx_raw, &R.cpx_chk}) if (o->a >= 0) for (auto& b : Bx) if (b.a == o->a && b.e == o->e) { o->se = se_cpx(b); break; }
  }
  if (cf.reserved_band == 1) {  
    auto on_band = [&](const Est& src, Est& dst, bool checked) {
      dst = Est(); dst.mask = src.mask; dst.status = src.status;
      if (src.a < 0 || src.status != OK) return;
      const Band* bb = nullptr; for (auto& b : Bc) if (b.a == src.a && b.e == src.e) { bb = &b; break; }
      if (!bb || !bb->valid) { dst.status = REF_DENOMINATOR; return; }
      dst.a = bb->a; dst.e = bb->e; dst.p = src.p; dst.rank = src.rank; dst.den = bb->den; dst.est = bb->est;
      if (checked) dst.status = !positive_denominator(*bb) ? REF_DENOMINATOR : !within_range(*bb) ? REF_RANGE : OK;
      if (dst.status != OK) dst.est = NaN;
      fin_cos(dst);
    };
    on_band(R.cpx_raw, R.cos_raw, false); on_band(R.cpx_chk, R.cos_chk, true);
    R.pjoint_on_cos = NaN;
    if (R.cos_chk.a >= 0) for (auto& b : Bc) if (b.a == R.cos_chk.a && b.e == R.cos_chk.e) { R.pjoint_on_cos = b.pj; break; }
  }
  
  if (std::isfinite(cf.v2)) {
    vector<double> cc(M), ssn(M);
    for (int k = 0; k < M; ++k) { const double t = (k + 1) * D; cc[k] = std::cos(t * cf.v2); ssn[k] = std::sin(t * cf.v2); }
    for (int cplx = 0; cplx < 2; ++cplx) {
      vector<Band> Br = B; Est& out = cplx ? R.rat_cpx : R.rat_cos; Est rawdummy;
      
      
      auto rows = [&](const double* c1, const double* s1, const double* cy, const double* sy, int k, double X[2][3], double y[2]) -> int {
        return rat_rows(cplx, cc.data(), ssn.data(), c1, s1, cy, sy, k, X, y);
      };
      auto fit = [&](const double* c1, const double* s1, const double* cy, const double* sy, int a, int e, double th[3]) {
        return rat_fit(cplx, cc.data(), ssn.data(), c1, s1, cy, sy, a, e, th);
      };
      
      auto coefs = [&](const double th[3], int k, double co[2][4]) -> int {
        const double mu = th[0], ka = th[1], rh = th[2];
        if (!cplx) { co[0][0] = 1 + rh * cc[k]; co[0][1] = 0; co[0][2] = -(mu + ka * cc[k]); co[0][3] = 0; return 1; }
        const double zr = cc[k], zi = ssn[k];
        co[0][0] = 1 + rh * zr; co[0][1] = -rh * zi; co[0][2] = -(mu + ka * zr); co[0][3] = ka * zi;
        co[1][0] = rh * zi; co[1][1] = 1 + rh * zr; co[1][2] = -ka * zi; co[1][3] = -(mu + ka * zr);
        return 2;
      };
      const int cidx4[4] = {iCY, iSY, iC1, iS1};
      vector<double> sdv((size_t)2 * M);
      vector<double> c1s(M), s1s(M), cys(M), sys(M);
      vector<double> thetas((size_t)Br.size() * 3, NaN);
      for (size_t bi = 0; bi < Br.size(); ++bi) {
        Band& b = Br[bi]; if (b.m < 6) continue;
        double th[3]; if (!fit(cosine_mean.data(), S1.data(), CY.data(), SY.data(), b.a, b.e, th)) continue;
        b.valid = true; b.mu = th[0]; b.est = th[0]; b.den = 1; for (int i = 0; i < 3; ++i) thetas[bi * 3 + i] = th[i];
        double T = 0;
        for (int k = b.a; k <= b.e; ++k) {
          double co[2][4]; const int nr = coefs(th, k, co);
          const double obs[4] = {CY[k], SY[k], cosine_mean[k], S1[k]};
          for (int q = 0; q < nr; ++q) {
            double v = 0, r = 0;
            for (int i = 0; i < 4; ++i) { r += co[q][i] * obs[i]; for (int j = 0; j < 4; ++j) if (co[q][i] != 0 && co[q][j] != 0) v += co[q][i] * co[q][j] * VB(cidx4[i], k, cidx4[j], k); }
            sdv[(size_t)q * M + k] = std::sqrt(std::max(v, 1e-300)); T = std::max(T, std::fabs(r) / sdv[(size_t)q * M + k]);
          }
        }
        int ex = 0;
        if (drawB) for (int r = 0; r < Rb; ++r) {
          for (int k = b.a; k <= b.e; ++k) {  
            const double zr = cplx ? cc[k] : cc[k], zi = cplx ? ssn[k] : 0.0;
            const double nr_ = th[0] + th[1] * zr, ni = th[1] * zi, dr = 1 + th[2] * zr, di = th[2] * zi, dd = dr * dr + di * di;
            const double qr = (nr_ * dr + ni * di) / dd, qi = (ni * dr - nr_ * di) / dd;  
            c1s[k] = cosine_mean[k] + zB(r, iC1, k); s1s[k] = S1[k] + zB(r, iS1, k);
            if (!cplx) { cys[k] = cosine_mean[k] * (th[0] + th[1] * cc[k]) / (1 + th[2] * cc[k]) + zB(r, iCY, k); sys[k] = SY[k]; }
            else { cys[k] = cosine_mean[k] * qr - S1[k] * qi + zB(r, iCY, k); sys[k] = cosine_mean[k] * qi + S1[k] * qr + zB(r, iSY, k); }
          }
          double ts[3]; if (!fit(c1s.data(), s1s.data(), cys.data(), sys.data(), b.a, b.e, ts)) continue;
          for (int k = b.a; k <= b.e; ++k) {
            double co[2][4]; const int nr = coefs(ts, k, co); const double obs[4] = {cys[k], sys[k], c1s[k], s1s[k]}; bool hit = false;
            for (int q = 0; q < nr; ++q) { double rr = 0; for (int i = 0; i < 4; ++i) rr += co[q][i] * obs[i]; if (std::fabs(rr) / sdv[(size_t)q * M + k] >= T) { hit = true; break; } }
            if (hit) { ++ex; break; }
          }
        }
        b.p = drawB ? (1.0 + ex) / (Rb + 1.0) : NaN;
      }
      auto se_rat = [&](const Band& b) {
        size_t bi = 0; for (; bi < Br.size(); ++bi) if (Br[bi].a == b.a && Br[bi].e == b.e) break;
        const double* th = &thetas[bi * 3];
        double A[3][3] = {{0}}; vector<std::array<double, 3>> Xr; vector<std::array<double, 4>> Cr; vector<int> kr;
        for (int k = b.a; k <= b.e; ++k) {
          double X[2][3], y[2]; const int nr = rows(cosine_mean.data(), S1.data(), CY.data(), SY.data(), k, X, y); double co[2][4]; coefs(th, k, co);
          for (int q = 0; q < nr; ++q) { Xr.push_back({X[q][0], X[q][1], X[q][2]}); Cr.push_back({co[q][0], co[q][1], co[q][2], co[q][3]}); kr.push_back(k);
            for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) A[i][j] += X[q][i] * X[q][j]; }
        }
        double Ai[3][3]; if (!inv3(A, Ai)) return NaN;
        double Mid[3][3] = {{0}};
        for (size_t u = 0; u < Xr.size(); ++u) for (size_t w = 0; w < Xr.size(); ++w) {
          double om_ = 0; for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) if (Cr[u][i] != 0 && Cr[w][j] != 0) om_ += Cr[u][i] * Cr[w][j] * VB(cidx4[i], kr[u], cidx4[j], kr[w]);
          for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) Mid[i][j] += Xr[u][i] * Xr[w][j] * om_;
        }
        double v = 0; for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) v += Ai[0][i] * Mid[i][j] * Ai[j][0];
        return std::sqrt(std::max(v, 0.0));
      };
      auto positive_rational_denominator = [&](const Band& b) { return b.valid; };
      select_ranked_windows(Br, [](const Band& b) { return b.p; }, positive_rational_denominator, within_range, cf.pthr, R.frequency_band_visible, rawdummy, out);
      if (out.a >= 0) for (auto& b : Br) if (b.a == out.a && b.e == out.e) { out.se = se_rat(b); break; }
    }
  }
  return R;
}
}  
