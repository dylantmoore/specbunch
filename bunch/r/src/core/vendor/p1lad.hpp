// p1lad.hpp: verbatim copy of cpp/part1/src/part1_ll.cpp lines 30-102 (LAD line by Wesolowsky pivot descent, OLS start,
// fit_lad), for the PUF library p1lib.cpp (2026-11-06). Needs llcore.hpp (NG, BF) included first. Do not edit: edit part1_ll.cpp
// and re-copy.
#pragma once
#include <utility>
// ---------------------------------------------------------------- LAD line
// min_{a,b} sum_j |y_j - a - b x_j| by Wesolowsky's (1981) pivot descent: for a pivot point p the best slope
// through it is the weighted median of (y_j - y_p)/(x_j - x_p) with weights |x_j - x_p|; the new line passes
// through p and the median point, which becomes a pivot candidate. Stops when no point on the current line
// gives a strict decrease. Start: the point with the smallest |residual| from (a0, b0).
static inline double lad_obj(const double* x, const double* y, int m, double a, double b) {
  double f = 0;
  for (int j = 0; j < m; ++j) f += std::fabs(y[j] - a - b * x[j]);
  return f;
}
static inline bool lad_through(const double* x, const double* y, int m, int p, double& a, double& b, int& q) {
  std::pair<double, int> sl[NG];
  double wv[NG];
  int cnt = 0; double W = 0;
  for (int j = 0; j < m; ++j) {
    double dx = x[j] - x[p];
    if (dx == 0.0) continue;
    sl[cnt] = {(y[j] - y[p]) / dx, j};
    ++cnt;
  }
  if (cnt == 0) return false;
  std::sort(sl, sl + cnt);
  for (int i = 0; i < cnt; ++i) { wv[i] = std::fabs(x[sl[i].second] - x[p]); W += wv[i]; }
  double cum = 0; int pick = cnt - 1;
  for (int i = 0; i < cnt; ++i) { cum += wv[i]; if (cum >= 0.5 * W) { pick = i; break; } }
  b = sl[pick].first; q = sl[pick].second; a = y[p] - b * x[p];
  return std::isfinite(a) && std::isfinite(b);
}
static bool lad_line(const double* x, const double* y, int m, double a0, double b0, double& a, double& b) {
  if (m < 2) return false;
  int p = 0; double best = INFINITY;
  for (int j = 0; j < m; ++j) { double r = std::fabs(y[j] - a0 - b0 * x[j]); if (r < best) { best = r; p = j; } }
  int q;
  if (!lad_through(x, y, m, p, a, b, q)) return false;
  double f = lad_obj(x, y, m, a, b);
  for (int it = 0; it < 500; ++it) {
    bool improved = false;
    for (int l = 0; l < m && !improved; ++l) {
      if (!(std::fabs(y[l] - a - b * x[l]) <= 1e-10 * (1.0 + std::fabs(y[l])))) continue;  // points on the line
      double a2, b2; int q2;
      if (!lad_through(x, y, m, l, a2, b2, q2)) continue;
      double f2 = lad_obj(x, y, m, a2, b2);
      if (f2 < f - 1e-13 * (1.0 + f)) { a = a2; b = b2; f = f2; improved = true; }
    }
    if (!improved) return true;
  }
  return true;
}

// OLS start for LAD (unweighted), on the given arrays
static inline void ols_line(const double* x, const double* y, int m, double& a, double& b) {
  double sx = 0, sy = 0;
  for (int j = 0; j < m; ++j) { sx += x[j]; sy += y[j]; }
  double xb = sx / m, yb = sy / m, Sxx = 0, Sxy = 0;
  for (int j = 0; j < m; ++j) { double dx = x[j] - xb; Sxx += dx * dx; Sxy += dx * (y[j] - yb); }
  b = Sxx > 0 ? Sxy / Sxx : 0.0; a = yb - b * xb;
}

// LAD band fit (same validity rules as F1): invalid if any C_j <= 0 or slope <= 0.
// warm: if finite, (a0, b0) start; else OLS start.
static inline bool fit_lad(const double* C, const double* x, int m, double a0, double b0, BF& r) {
  double y[NG];
  for (int j = 0; j < m; ++j) {
    if (!(C[j] > 0.0)) return false;
    y[j] = std::log(C[j]);
  }
  if (!(std::isfinite(a0) && std::isfinite(b0))) ols_line(x, y, m, a0, b0);
  double a, b;
  if (!lad_line(x, y, m, a0, b0, a, b)) return false;
  if (!(b > 0.0) || !std::isfinite(a) || !std::isfinite(b)) return false;
  r.logB = a; r.B = std::exp(a); r.sigma = std::sqrt(b); r.logS = 0.5 * std::log(b);
  return std::isfinite(r.B) && std::isfinite(r.sigma);
}
