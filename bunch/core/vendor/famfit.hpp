// famfit.hpp: levels Gaussian-QML (weighted least squares) fits of the cosine ECF for the Cauchy and Voigt mean
// families, the generic sup-residual statistic, and a runtime-width copy of llcore's multiplier-bootstrap sums.
//
//   G: mu(t) = B exp(-sigma^2 t^2 / 2)            (llcore.hpp gqml_nr / fit_lev, used unchanged for the G family)
//   C: mu(t) = B exp(-a t)                        (cfit: gqml_nr transcribed to (log B, log a))
//   V: mu(t) = B exp(-a t - s2 t^2 / 2), a, s2 >= 0 (vfit: Newton in (log B, a', s2'), boundary fits on exit)
//
// All fits maximise L = -1/2 sum_j w_j (C_j - mu_j)^2 with w_j = n / Var_j (the F3 objective).
// Scale summary for every family: sigma_eq = 1 / t_e, where t_e solves a t + s2 t^2 / 2 = 1/2 (the frequency at
// which the friction CF falls to e^{-1/2}): sigma_eq = a + sqrt(a^2 + s2). G: sigma_eq = sigma. C: sigma_eq = 2a.
#pragma once
#include "llcore.hpp"

constexpr int GMAX = 160;  // max grid points for the C/V (extended) grid

static inline double wgt(double v, double nd) { return v > 1e-12 ? nd / v : nd; }

// ---------------------------------------------------------------- G on any grid width (gqml_nr, dynamic arrays)
// Identical arithmetic to llcore gqml_nr; used only as the a = 0 boundary of V on the extended grid.
static NR gfit_dyn(const double* C, const double* t, const double* var, int M, size_t n, double sigma_init,
                   int max_iter = 100, double tol = 1e-8) {
  double w[GMAX];
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

// ---------------------------------------------------------------- C: gqml_nr transcribed to mu = exp(alpha - a t)
// theta = (alpha, psi = log a), b_j = a t_j. Same full-Hessian / Gauss-Newton fallback, same step and stopping rule.
// Start: OLS of log max(C, 1e-10) on (1, -t); a_start = OLS slope if positive and finite, else a_init.
struct NRc { double B, a; bool conv; };
static NRc cfit(const double* C, const double* t, const double* var, int M, size_t n, double a_init,
                int max_iter = 100, double tol = 1e-8) {
  double w[GMAX];
  for (int j = 0; j < M; ++j) w[j] = var[j] > 1e-12 ? (double)n / var[j] : (double)n;
  double sxx = 0, sx = 0, sy = 0, sxy = 0;
  for (int j = 0; j < M; ++j) {
    double x = -t[j], y = std::log(std::max(C[j], 1e-10));
    sx += x; sy += y; sxx += x * x; sxy += x * y;
  }
  double Md = (double)M, denom = Md * sxx - sx * sx, alpha, ai;
  if (std::fabs(denom) > 1e-30) {
    double b1 = (Md * sxy - sx * sy) / denom;
    alpha = (sy - b1 * sx) / Md;
    ai = (b1 > 0.0 && std::isfinite(b1)) ? b1 : a_init;
  } else { alpha = std::log(std::max(C[0], 1e-10)); ai = a_init; }
  double psi = std::log(ai);
  bool conv = false;
  for (int it = 0; it < max_iter; ++it) {
    double a = std::exp(psi);
    double ga = 0, gp = 0, Haa = 0, Hap = 0, Hpp = 0, Ga = 0, Gap = 0, Gpp = 0;
    for (int j = 0; j < M; ++j) {
      double bj = a * t[j];
      double eta = std::clamp(alpha - bj, -500.0, 500.0), mj = std::exp(eta);
      double rj = C[j] - mj, wj = w[j], wm = wj * mj, wm2 = wm * mj;
      ga += wm * rj; gp -= wm * bj * rj;
      Haa += wm * (rj - mj); Hap += wm * bj * (mj - rj);
      Hpp += wm * bj * ((bj - 1.0) * rj - bj * mj);
      Ga -= wm2; Gap += wm2 * bj; Gpp -= wm2 * bj * bj;
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

// ---------------------------------------------------------------- V: constrained (a, s2 >= 0)
// Interior: Newton in theta = (alpha, a', s2'), a' = a t_s, s2' = s2 t_s^2, t_s = band midpoint (so a', s2' are O(1)).
// Full Hessian when -H is positive definite, else Gauss-Newton; step halving until L does not decrease.
// Kinds: 0 interior, 1 boundary a = 0 (G), 2 boundary s2 = 0 (C), 3 interior reached by the KKT restart.
// Rule: start at the OLS fit of log max(C, 1e-10) on (1, -t, -t^2/2) if it is feasible. If that start is infeasible,
// an iterate leaves the feasible set, or Newton fails to converge, fit both boundary models and keep the one with the
// larger L. KKT: if L still increases into the interior from that boundary point (dL/da > 0 at a = 0, or dL/ds2 > 0
// at s2 = 0), restart the interior Newton there (moved 5% of the scale inward) and keep it if it converges feasibly
// with a larger L. With the boundary models solved to their own optimum, this returns a KKT point of the
// constrained problem.
struct VF { double B = NaN, a = NaN, s2 = NaN, L = -INFINITY; int kind = -1; bool ok = false; };

static inline double Lval(const double* C, const double* t, const double* w, int M, double al, double a, double s2) {
  double L = 0;
  for (int j = 0; j < M; ++j) {
    double eta = std::clamp(al - a * t[j] - 0.5 * s2 * t[j] * t[j], -500.0, 500.0), r = C[j] - std::exp(eta);
    L -= 0.5 * w[j] * r * r;
  }
  return L;
}

static bool solve3(const double A[3][3], const double b[3], double x[3], bool& pd) {
  // Cholesky; pd = false if not PD.
  double L[3][3] = {{0}};
  pd = true;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j <= i; ++j) {
      double s = A[i][j];
      for (int k = 0; k < j; ++k) s -= L[i][k] * L[j][k];
      if (i == j) { if (!(s > 1e-300)) { pd = false; return false; } L[i][i] = std::sqrt(s); }
      else L[i][j] = s / L[j][j];
    }
  double y[3];
  for (int i = 0; i < 3; ++i) { double s = b[i]; for (int k = 0; k < i; ++k) s -= L[i][k] * y[k]; y[i] = s / L[i][i]; }
  for (int i = 2; i >= 0; --i) { double s = y[i]; for (int k = i + 1; k < 3; ++k) s -= L[k][i] * x[k]; x[i] = s / L[i][i]; }
  return std::isfinite(x[0]) && std::isfinite(x[1]) && std::isfinite(x[2]);
}

// returns 0 converged feasible, 1 left the feasible set, 2 failed. tau = t / t_s.
static int newton3(const double* C, const double* tau, const double* w, int M, double th[3],
                   int max_iter = 100, double tol = 1e-8) {
  for (int it = 0; it < max_iter; ++it) {
    double g[3] = {0, 0, 0}, H[3][3] = {{0}}, G[3][3] = {{0}};
    for (int j = 0; j < M; ++j) {
      double x0 = 1.0, x1 = -tau[j], x2 = -0.5 * tau[j] * tau[j];
      double eta = std::clamp(th[0] + th[1] * x1 + th[2] * x2, -500.0, 500.0), m = std::exp(eta);
      double r = C[j] - m, wm = w[j] * m, c = wm * (m - r), d = wm * m, x[3] = {x0, x1, x2};
      for (int p = 0; p < 3; ++p) {
        g[p] += wm * r * x[p];
        for (int q = 0; q <= p; ++q) { H[p][q] += c * x[p] * x[q]; G[p][q] += d * x[p] * x[q]; }
      }
    }
    for (int p = 0; p < 3; ++p) for (int q = p + 1; q < 3; ++q) { H[p][q] = H[q][p]; G[p][q] = G[q][p]; }
    double dlt[3]; bool pd;
    if (!solve3(H, g, dlt, pd)) { if (!solve3(G, g, dlt, pd)) return 2; }
    double L0 = Lval(C, tau, w, M, th[0], th[1], th[2]);
    double st = 1.0, nt[3];
    for (int h = 0; h < 40; ++h) {
      for (int p = 0; p < 3; ++p) nt[p] = th[p] + st * dlt[p];
      double L1 = Lval(C, tau, w, M, nt[0], nt[1], nt[2]);
      if (L1 >= L0 || !std::isfinite(L0)) break;
      st *= 0.5;
    }
    for (int p = 0; p < 3; ++p) th[p] = nt[p];
    if (!(std::isfinite(th[0]) && std::isfinite(th[1]) && std::isfinite(th[2]))) return 2;
    if (th[1] < 0.0 || th[2] < 0.0) return 1;
    double dd = dlt[0] * dlt[0] + dlt[1] * dlt[1] + dlt[2] * dlt[2];
    if (dd < tol * tol * (1.0 + th[0] * th[0] + th[1] * th[1] + th[2] * th[2])) return 0;
  }
  return 2;
}

struct VDiag { int restarts = 0; };
static VF vfit(const double* C, const double* t, const double* var, int M, size_t n, double sig_init, VDiag* dg = nullptr) {
  double w[GMAX], tau[GMAX];
  const double ts = 0.5 * (t[0] + t[M - 1]);
  for (int j = 0; j < M; ++j) { w[j] = wgt(var[j], (double)n); tau[j] = t[j] / ts; }
  VF best;
  // OLS start
  double A[3][3] = {{0}}, b[3] = {0, 0, 0};
  for (int j = 0; j < M; ++j) {
    double x[3] = {1.0, -tau[j], -0.5 * tau[j] * tau[j]}, y = std::log(std::max(C[j], 1e-10));
    for (int p = 0; p < 3; ++p) { b[p] += x[p] * y; for (int q = 0; q < 3; ++q) A[p][q] += x[p] * x[q]; }
  }
  double th[3]; bool pd;
  bool st_ok = solve3(A, b, th, pd) && th[1] >= 0.0 && th[2] >= 0.0;
  if (st_ok && newton3(C, tau, w, M, th) == 0) {
    best.B = std::exp(std::clamp(th[0], -500.0, 500.0)); best.a = th[1] / ts; best.s2 = th[2] / (ts * ts);
    best.L = Lval(C, tau, w, M, th[0], th[1], th[2]); best.kind = 0; best.ok = true;
    return best;
  }
  // boundary fits
  NR g = gfit_dyn(C, t, var, M, n, sig_init);
  NRc c = cfit(C, t, var, M, n, 0.5 * sig_init);
  VF fg, fc;
  if (g.conv && std::isfinite(g.B) && std::isfinite(g.sigma)) {
    fg.B = g.B; fg.a = 0.0; fg.s2 = g.sigma * g.sigma; fg.kind = 1; fg.ok = true;
    fg.L = Lval(C, t, w, M, std::log(g.B), 0.0, fg.s2);
  }
  if (c.conv && std::isfinite(c.B) && std::isfinite(c.a)) {
    fc.B = c.B; fc.a = c.a; fc.s2 = 0.0; fc.kind = 2; fc.ok = true;
    fc.L = Lval(C, t, w, M, std::log(c.B), c.a, 0.0);
  }
  if (fg.ok && (!fc.ok || fg.L >= fc.L)) best = fg; else if (fc.ok) best = fc; else return best;
  // KKT check at the chosen boundary point (scaled coordinates)
  {
    const double al = std::log(best.B), ap = best.a * ts, sp = best.s2 * ts * ts;
    double ga = 0, gs = 0;
    for (int j = 0; j < M; ++j) {
      double eta = std::clamp(al - ap * tau[j] - 0.5 * sp * tau[j] * tau[j], -500.0, 500.0), m = std::exp(eta);
      double wr = w[j] * (C[j] - m) * m;
      ga -= wr * tau[j]; gs -= wr * 0.5 * tau[j] * tau[j];
    }
    bool viol = (best.kind == 1) ? (ga > 0.0) : (gs > 0.0);
    if (viol) {
      if (dg) dg->restarts++;
      double t2[3] = {al, ap, sp};
      if (best.kind == 1) t2[1] = 0.05 * std::sqrt(sp); else t2[2] = (0.05 * ap) * (0.05 * ap);
      if (newton3(C, tau, w, M, t2) == 0) {
        double L2 = Lval(C, tau, w, M, t2[0], t2[1], t2[2]);
        if (L2 > best.L) {
          best.B = std::exp(std::clamp(t2[0], -500.0, 500.0)); best.a = t2[1] / ts; best.s2 = t2[2] / (ts * ts);
          best.L = L2; best.kind = 3;
        }
      }
    }
  }
  return best;
}

// ---------------------------------------------------------------- family band fit
// fam: 1 = C, 2 = V. Output in the BF struct (sigma = sigma_eq, logS = log sigma_eq, seS = se(log sigma_eq)),
// plus a, s2, kind. SE: inverse Gauss-Newton weighted information (as fit_lev), delta method for log sigma_eq.
struct BFX { BF r; double a = NaN, s2 = NaN; int kind = -1; };
static inline double sig_eq(double a, double s2) { return a + std::sqrt(a * a + s2); }
static inline double mu_v(double B, double a, double s2, double t) { return B * std::exp(-a * t - 0.5 * s2 * t * t); }

template <bool WANT_SE>
static bool fit_fam(int fam, const double* C, const double* t, const double* var, int m, size_t n, double sig_init,
                    BFX& o, VDiag* dg = nullptr) {
  BF& r = o.r;
  if (fam == 1) {
    NRc q = cfit(C, t, var, m, n, 0.5 * sig_init);
    if (!(q.conv && std::isfinite(q.B) && std::isfinite(q.a))) return false;
    o.a = q.a; o.s2 = 0.0; o.kind = 2;
    r.B = q.B; r.sigma = 2.0 * q.a; r.logB = std::log(q.B); r.logS = std::log(r.sigma);
    if (WANT_SE) {
      double J11 = 0, J12 = 0, J22 = 0;
      for (int j = 0; j < m; ++j) {
        double wj = wgt(var[j], (double)n), bj = q.a * t[j], mj = q.B * std::exp(-bj), wm2 = wj * mj * mj;
        J11 += wm2; J12 -= wm2 * bj; J22 += wm2 * bj * bj;
      }
      double det = J11 * J22 - J12 * J12;
      r.seB = std::sqrt(J22 / det); r.seS = std::sqrt(J11 / det);
      r.sevalid = det > 0 && std::isfinite(r.seB) && std::isfinite(r.seS) && r.seB > 0 && r.seS > 0;
    }
  } else {
    VF v = vfit(C, t, var, m, n, sig_init, dg);
    if (!(v.ok && std::isfinite(v.B) && v.B > 0 && std::isfinite(v.a) && std::isfinite(v.s2))) return false;
    o.a = v.a; o.s2 = v.s2; o.kind = v.kind;
    r.B = v.B; r.sigma = sig_eq(v.a, v.s2); r.logB = std::log(v.B); r.logS = std::log(r.sigma);
    if (!(std::isfinite(r.sigma) && r.sigma > 0)) return false;
    if (WANT_SE) {
      if (v.kind == 1 || v.kind == 2) {
        // boundary: the boundary model's own information (G in (log B, log sigma), C in (log B, log a))
        double J11 = 0, J12 = 0, J22 = 0;
        for (int j = 0; j < m; ++j) {
          double wj = wgt(var[j], (double)n);
          double bj = (v.kind == 1) ? v.s2 * t[j] * t[j] : v.a * t[j];
          double mj = mu_v(v.B, v.a, v.s2, t[j]), wm2 = wj * mj * mj;
          J11 += wm2; J12 -= wm2 * bj; J22 += wm2 * bj * bj;
        }
        double det = J11 * J22 - J12 * J12;
        r.seB = std::sqrt(J22 / det); r.seS = std::sqrt(J11 / det);
        r.sevalid = det > 0 && std::isfinite(r.seB) && std::isfinite(r.seS) && r.seB > 0 && r.seS > 0;
      } else {
        const double ts = 0.5 * (t[0] + t[m - 1]);
        double J[3][3] = {{0}};
        for (int j = 0; j < m; ++j) {
          double tau = t[j] / ts, x[3] = {1.0, -tau, -0.5 * tau * tau};
          double mj = mu_v(v.B, v.a, v.s2, t[j]), wm2 = wgt(var[j], (double)n) * mj * mj;
          for (int p = 0; p < 3; ++p) for (int q = 0; q < 3; ++q) J[p][q] += wm2 * x[p] * x[q];
        }
        // invert J (3x3) via cofactors
        double c00 = J[1][1] * J[2][2] - J[1][2] * J[2][1], c01 = -(J[1][0] * J[2][2] - J[1][2] * J[2][0]),
               c02 = J[1][0] * J[2][1] - J[1][1] * J[2][0];
        double det = J[0][0] * c00 + J[0][1] * c01 + J[0][2] * c02;
        double I[3][3];
        I[0][0] = c00 / det; I[0][1] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) / det; I[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / det;
        I[1][0] = c01 / det; I[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / det; I[1][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) / det;
        I[2][0] = c02 / det; I[2][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) / det; I[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / det;
        double ap = v.a * ts, sp = v.s2 * ts * ts, R = std::sqrt(ap * ap + sp), se = ap + R;
        double g1 = (1.0 + ap / R) / se, g2 = 0.5 / (R * se);
        double vS = g1 * g1 * I[1][1] + 2 * g1 * g2 * I[1][2] + g2 * g2 * I[2][2];
        r.seB = std::sqrt(I[0][0]); r.seS = std::sqrt(vS);
        r.sevalid = det > 0 && std::isfinite(r.seB) && std::isfinite(r.seS) && r.seB > 0 && r.seS > 0;
      }
    }
  }
  return std::isfinite(r.B) && std::isfinite(r.sigma);
}

static double sup_std_resid_v(const double* C, const double* t, const double* s, int M, double B, double a, double s2) {
  double D = 0;
  for (int j = 0; j < M; ++j) { double d = std::fabs(C[j] - mu_v(B, a, s2, t[j])) / s[j]; if (d > D) D = d; }
  return D;
}

// ---------------------------------------------------------------- multiplier-bootstrap sums, runtime width
// Copy of llcore boot_sums with M columns padded to MPd (multiple of 16): same signs, chunking and accumulation
// order. Columns j < 77 agree with llcore boot_sums to ~1e-12 but NOT bitwise (the compiled per-observation
// recurrence differs), so the G family always uses llcore boot_sums itself.
static void boot_sums_w(const vector<double>& dZ, double t0, double Delta, uint64_t seed, int M, int MPd,
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
#pragma omp for schedule(dynamic, 1)
    for (int ch = 0; ch < NCH; ++ch) {
      double* Yl = &Ych[(size_t)ch * YS]; double* Sl = &Sch[(size_t)ch * BOOT_B];
      const size_t g0 = ch * per, g1 = std::min(ng, g0 + per);
      for (size_t g = g0; g < g1; ++g) {
        size_t i0 = g * 64; int L = (int)std::min<size_t>(64, n - i0);
        for (int ii = 0; ii < 64; ++ii) {
          double* X = &Xv[(size_t)ii * MPd];
          if (ii >= L) { for (int j = 0; j < MPd; ++j) X[j] = 0.0; continue; }
          double z = dZ[i0 + ii];
          double cl = std::cos(t0 * z), sl = std::sin(t0 * z), cd = std::cos(Delta * z), sd = std::sin(Delta * z);
          X[0] = cl; double cp = cl, sp = sl;
          for (int j = 1; j < M; ++j) {
            double cn = cp * cd - sp * sd, sn = sp * cd + cp * sd;
            X[j] = cn; cp = cn; sp = sn;
          }
          for (int j = M; j < MPd; ++j) X[j] = 0.0;
        }
        const uint64_t mask = (L == 64) ? ~0ULL : ((1ULL << L) - 1);
        for (int b = 0; b < BOOT_B; ++b) {
          uint64_t w = sign_word(key, b, g);
          Sl[b] += 2.0 * __builtin_popcountll(w & mask) - (double)L;
          double* yb = &Yl[(size_t)b * MPd];
          for (int jc = 0; jc < MPd; jc += 16) {
            double acc[16];
            for (int q = 0; q < 16; ++q) acc[q] = yb[jc + q];
            for (int ii = 0; ii < 64; ++ii) {
              double s = ((w >> ii) & 1ULL) ? 1.0 : -1.0;
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
