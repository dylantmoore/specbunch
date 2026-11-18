








#pragma once
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace sbpol {

struct Cand {
  int rank, a, e;
  double B, sigma, seB, seS, p;
  bool fit_converged, coherent_values, fit_test_passed;
  std::string win;
};


inline Cand make_cand(int rank, int a, int e, double B, double s, double seB, double seS, double p, double rok) {
  Cand c;
  c.rank = rank; c.a = a; c.e = e; c.B = B; c.sigma = s; c.seB = seB; c.seS = seS; c.p = p;
  const bool fin = std::isfinite(B) && std::isfinite(s);
  c.fit_converged = fin;
  c.coherent_values = fin && 0 < B && B < 1 && s > 0;
  c.fit_test_passed = std::isfinite(p) && p >= 0.10;
  char buf[64];
  std::snprintf(buf, sizeof buf, "u%.2f-%.2f", 0.20 + 0.05 * a, 0.20 + 0.05 * e);
  c.win = buf;
  return c;
}

inline std::string rle(const std::vector<std::string>& xs, size_t upto) {
  std::string out; std::string prev; size_t k = 0; bool have = false;
  auto flush = [&]() {
    if (!have) return;
    if (!out.empty()) out += ';';
    out += k == 1 ? prev : prev + "*" + std::to_string(k);
  };
  for (size_t i = 0; i < upto && i < xs.size(); ++i) {
    if (have && xs[i] == prev) { ++k; continue; }
    flush(); prev = xs[i]; k = 1; have = true;
  }
  flush();
  return out;
}


inline int depth_index(const std::string& s) {
  static const char* D[9] = {"fit_converged", "coherent_values", "cdf_inconsistent", "moments_inconsistent", "no_valid_distribution", "residual_cdf_inconsistent", "distribution_check_inconclusive", "poor_fit", "selected"};
  for (int i = 0; i < 9; ++i) if (s == D[i]) return i;
  return -1;
}

struct Stage {
  bool answered = false;
  int idx = -1;
  double B = NAN, sigma = NAN;
  int rank = -1;
  std::string win, reason, higher;
};

struct PolicyOut {
  Stage st[4];                      
  std::vector<std::string> labels;  
  int n_screened = 0;
};



inline PolicyOut policy(const std::vector<Cand>& cands, const std::string& internal,
                        const std::function<std::string(int)>& distribution_checks) {
  PolicyOut R;
  const int nc = (int)cands.size();
  if (!internal.empty() || nc == 0) {
    for (auto& s : R.st) { s.answered = false; s.reason = "no_coherent_window"; s.higher = ""; }
    return R;
  }
  std::vector<int> memo_has(nc, 0);
  std::vector<std::string> memo(nc);
  auto check_distribution = [&](int i) -> const std::string& {
    if (!memo_has[i]) { memo[i] = distribution_checks ? distribution_checks(i) : std::string(); memo_has[i] = 1; }
    return memo[i];
  };
  auto need = [&](int k, const Cand& c) {
    if (!c.fit_converged) return false;
    if (k >= 1 && !c.coherent_values) return false;
    if (k >= 3 && !c.fit_test_passed) return false;
    return true;
  };
  int stop[4];
  for (int k = 0; k < 4; ++k) {
    Stage& s = R.st[k];
    stop[k] = nc;
    for (int i = 0; i < nc; ++i) {
      if (!need(k, cands[i])) continue;
      if (k < 2) { s.answered = true; s.idx = i; break; }
      if (check_distribution(i).empty()) { s.answered = true; s.idx = i; break; }
    }
  }
  
  for (int k = 0; k < 4; ++k) {
    Stage& s = R.st[k];
    if (s.answered) continue;
    if (k == 0) s.reason = "no_coherent_window";
    else if (k == 1) s.reason = R.st[0].answered ? "no_coherent_window" : "no_coherent_window";
    else if (k == 2) {
      int best = -1; std::string lab;
      for (int i = 0; i < nc; ++i)
        if (cands[i].fit_converged && cands[i].coherent_values && memo_has[i]) {
          const int d = depth_index(memo[i]);  
          if (d > best) { best = d; lab = memo[i]; }
        }
      s.reason = best >= 0 ? "no_window_passes:" + lab : R.st[1].reason;
    } else {
      const Stage& prev = R.st[k - 1];
      s.reason = prev.answered ? "no_window_passes:poor_fit" : prev.reason;
    }
  }
  
  R.labels.resize(nc);
  for (int i = 0; i < nc; ++i) {
    const Cand& c = cands[i];
    std::string L;
    if (!c.fit_converged) L = "fit_converged";
    else if (!c.coherent_values) L = "coherent_values";
    else if (memo_has[i]) L = !memo[i].empty() ? memo[i] : (!c.fit_test_passed ? "poor_fit" : "pass");
    else L = !c.fit_test_passed ? "poor_fit_distribution_checks_not_evaluated" : "distribution_checks_not_evaluated";
    R.labels[i] = L;
  }
  for (int i = 0; i < nc; ++i) R.n_screened += memo_has[i];
  for (int k = 0; k < 4; ++k) {
    Stage& s = R.st[k];
    if (s.answered) {
      const Cand& c = cands[s.idx];
      s.B = k >= 1 ? c.B : std::fmin(std::fmax(c.B, 0.0), 1.0);  
      s.sigma = c.sigma; s.rank = c.rank; s.win = c.win;
      s.higher = rle(R.labels, (size_t)s.idx);
    } else {
      s.higher = rle(R.labels, (size_t)stop[k]);
    }
  }
  return R;
}

}  
