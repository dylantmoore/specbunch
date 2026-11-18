















#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "Highs.h"

extern "C" {
int moment_box(int K, const double* c, const double* s, int n, double B, double sigma, double* lo, double* hi, double* diskmargin);
int cdf_screen(const double* x, const double* f, int n, double B, double sig, double* margin, int* which);
int probit_screen(const double* x, const double* f, int n, double B, double sigma, double* margin, int* endpoint);
}

namespace sbdistribution {

constexpr int K = 128, D = 256, J = 66;
constexpr double kPi = 3.141592653589793;  
using cplx = std::complex<double>;

struct Summ {
  int n = 0;
  double c[K], s[K], x[J], f[J];
};



inline void herm_eig(int N, std::vector<cplx> A, std::vector<double>& ev, std::vector<cplx>& V) {
  V.assign((size_t)N * N, cplx(0, 0));
  for (int i = 0; i < N; ++i) V[(size_t)i * N + i] = 1.0;
  auto a = [&](int i, int j) -> cplx& { return A[(size_t)i * N + j]; };
  for (int sweep = 0; sweep < 100; ++sweep) {
    double off = 0, diag = 0;
    for (int i = 0; i < N; ++i) {
      diag += std::norm(a(i, i));
      for (int j = i + 1; j < N; ++j) off += std::norm(a(i, j));
    }
    if (off <= 1e-32 * diag || off == 0.0) break;
    for (int p = 0; p < N - 1; ++p)
      for (int q = p + 1; q < N; ++q) {
        const cplx z = a(p, q);
        const double r = std::abs(z);
        if (r == 0.0) continue;
        const double app = a(p, p).real(), aqq = a(q, q).real();
        if (r < 1e-300 || (sweep > 3 && r <= 1e-18 * std::sqrt(std::fabs(app * aqq)))) { a(p, q) = 0; a(q, p) = 0; continue; }
        const cplx e = z / r, eb = std::conj(e);
        const double theta = (aqq - app) / (2.0 * r);
        const double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        const double cs = 1.0 / std::sqrt(t * t + 1.0), sn = t * cs;
        
        const cplx g00 = cs, g01 = sn, g10 = -sn * eb, g11 = cs * eb;
        for (int k = 0; k < N; ++k) {  
          const cplx akp = a(k, p), akq = a(k, q);
          a(k, p) = akp * g00 + akq * g10; a(k, q) = akp * g01 + akq * g11;
        }
        for (int k = 0; k < N; ++k) {  
          const cplx apk = a(p, k), aqk = a(q, k);
          a(p, k) = std::conj(g00) * apk + std::conj(g10) * aqk; a(q, k) = std::conj(g01) * apk + std::conj(g11) * aqk;
        }
        a(p, q) = 0; a(q, p) = 0; a(p, p) = a(p, p).real(); a(q, q) = a(q, q).real();
        for (int k = 0; k < N; ++k) {
          const cplx vkp = V[(size_t)k * N + p], vkq = V[(size_t)k * N + q];
          V[(size_t)k * N + p] = vkp * g00 + vkq * g10; V[(size_t)k * N + q] = vkp * g01 + vkq * g11;
        }
      }
  }
  std::vector<int> ord(N);
  for (int i = 0; i < N; ++i) ord[i] = i;
  std::stable_sort(ord.begin(), ord.end(), [&](int x, int y) { return a(x, x).real() < a(y, y).real(); });
  ev.resize(N);
  std::vector<cplx> W((size_t)N * N);
  for (int j = 0; j < N; ++j) {
    ev[j] = a(ord[j], ord[j]).real();
    double nn = 0;
    for (int k = 0; k < N; ++k) nn += std::norm(V[(size_t)k * N + ord[j]]);
    nn = std::sqrt(nn);
    for (int k = 0; k < N; ++k) W[(size_t)k * N + j] = V[(size_t)k * N + ord[j]] / nn;
  }
  V.swap(W);
}



inline std::string compat_screen(const Summ& S, double B, double sigma) {
  const double radius = std::sqrt(2 * std::log(512 / .01) / (double)S.n);
  double lo[2][K], hi[2][K], d[K];
  for (int h = 0; h < K; ++h) {
    lo[0][h] = std::max(-1.0, S.c[h] - radius); lo[1][h] = std::max(-1.0, S.s[h] - radius);
    hi[0][h] = std::min(1.0, S.c[h] + radius); hi[1][h] = std::min(1.0, S.s[h] + radius);
    const double fr = (double)(h + 1) / 4, v = sigma * fr;
    d[h] = std::exp(-.5 * (v * v));
  }
  double m = -INFINITY;
  for (int h = 0; h < K; ++h) {
    const double center = B * d[h], rad = (1 - B) * d[h];
    const double dx = std::max(std::max(lo[0][h] - center, center - hi[0][h]), 0.);
    const double dy = std::max(std::max(lo[1][h], -hi[1][h]), 0.);
    const double mg = std::hypot(dx, dy) - rad;
    if (mg > m) m = mg;
  }
  if (m > 1e-12) return "coordinate_exclusion";
  int Kp = 0;
  for (int h = 0; h < K; ++h) Kp += d[h] >= 1e-8;
  if (Kp == 0) return "unresolved";
  const int N = Kp + 1;
  std::vector<cplx> first(N);
  first[0] = cplx(1 - B, 0);
  for (int h = 0; h < Kp; ++h) first[h + 1] = cplx(lo[0][h] + hi[0][h], lo[1][h] + hi[1][h]) / (2 * d[h]) - B;
  std::vector<cplx> T((size_t)N * N);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) T[(size_t)i * N + j] = j >= i ? first[j - i] : std::conj(first[i - j]);
  std::vector<double> ev; std::vector<cplx> V;
  herm_eig(N, T, ev, V);
  const int count = std::min(5, Kp + 1);
  double best_norm = INFINITY, best_upper = 0, best_scale = 1;
  for (int jv = 0; jv < count; ++jv) {
    std::vector<cplx> v(N);
    for (int i = 0; i < N; ++i) v[i] = V[(size_t)i * N + jv];
    double norm = 0;
    for (int i = 0; i < N; ++i) norm += std::norm(v[i]);
    double sum_ar = 0, rt = 0, it = 0, sr = 0, si = 0;
    for (int h = 1; h <= Kp; ++h) {
      cplx acc(0, 0);
      for (int i = 0; i + h < N; ++i) acc += std::conj(v[i]) * v[i + h];
      const cplx ah = 2.0 * acc;
      sum_ar += ah.real();
      const double re = ah.real() / d[h - 1], im = -ah.imag() / d[h - 1];
      rt += std::max(re * lo[0][h - 1], re * hi[0][h - 1]);
      it += std::max(im * lo[1][h - 1], im * hi[1][h - 1]);
      sr += std::max(std::fabs(re * lo[0][h - 1]), std::fabs(re * hi[0][h - 1]));
      si += std::max(std::fabs(im * lo[1][h - 1]), std::fabs(im * hi[1][h - 1]));
    }
    const double constant = (1 - B) * norm - B * sum_ar;
    const double upper = constant + rt + it;
    const double scale = std::max(1., std::fabs(constant) + sr + si);
    if (upper / scale < best_norm) { best_norm = upper / scale; best_upper = upper; best_scale = scale; }
  }
  if (best_upper < -1e-9 * best_scale) return "joint_exclusion";
  return best_upper < 0 ? "unresolved" : "unexcluded";
}


struct LPRes { int status = 4; std::vector<double> x, ineq_marg; long long iters = 0; };



inline LPRes linprog_highs(int ncol, const std::vector<double>& cost, int nub, const std::vector<double>& Aub ,
                           const std::vector<double>& bub, const std::vector<double>& Aeq, double beq, long long iter_limit) {
  LPRes R;
  const int nrow = nub + 1;
  HighsLp lp;
  lp.num_col_ = ncol; lp.num_row_ = nrow;
  lp.col_cost_ = cost;
  lp.col_lower_.assign(ncol, 0.0); lp.col_upper_.assign(ncol, kHighsInf);
  lp.row_lower_.assign(nrow, -kHighsInf); lp.row_upper_.assign(nrow, 0.0);
  for (int i = 0; i < nub; ++i) lp.row_upper_[i] = bub[i];
  lp.row_lower_[nub] = beq; lp.row_upper_[nub] = beq;
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = ncol; lp.a_matrix_.num_row_ = nrow;
  lp.a_matrix_.start_.assign(ncol + 1, 0);
  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < nub; ++i) {
      const double v = Aub[(size_t)i * ncol + j];
      if (v != 0.0) { lp.a_matrix_.index_.push_back(i); lp.a_matrix_.value_.push_back(v); }
    }
    if (Aeq[j] != 0.0) { lp.a_matrix_.index_.push_back(nub); lp.a_matrix_.value_.push_back(Aeq[j]); }
    lp.a_matrix_.start_[j + 1] = (HighsInt)lp.a_matrix_.index_.size();
  }
  Highs h;
  h.setOptionValue("output_flag", false);
  h.setOptionValue("log_to_console", false);
  h.setOptionValue("presolve", "on");
  h.setOptionValue("highs_debug_level", 0);
  h.setOptionValue("dual_feasibility_tolerance", 1e-8);
  h.setOptionValue("primal_feasibility_tolerance", 1e-8);
  h.setOptionValue("simplex_strategy", 1);  
  h.setOptionValue("threads", 1);
  const HighsInt lim = (HighsInt)std::max<long long>(1, std::min<long long>(iter_limit, 2000000000LL));
  h.setOptionValue("simplex_iteration_limit", lim);
  h.setOptionValue("ipm_iteration_limit", lim);
  if (h.passModel(lp) == HighsStatus::kError) { R.status = 2; return R; }  
  const HighsStatus rs = h.run();
  const HighsModelStatus ms = h.getModelStatus();
  const HighsInfo& info = h.getInfo();
  R.iters = (long long)info.simplex_iteration_count + (long long)std::max<HighsInt>(0, info.ipm_iteration_count) +
            (long long)std::max<HighsInt>(0, info.crossover_iteration_count);
  auto scipy_status = [](HighsModelStatus s) {
    switch (s) {
      case HighsModelStatus::kOptimal: return 0;
      case HighsModelStatus::kTimeLimit: case HighsModelStatus::kIterationLimit: return 1;
      case HighsModelStatus::kModelError: case HighsModelStatus::kInfeasible: return 2;
      case HighsModelStatus::kUnbounded: return 3;
      default: return 4;
    }
  };
  if (rs == HighsStatus::kError || ms != HighsModelStatus::kOptimal) { R.status = scipy_status(ms); return R; }
  const HighsSolution& sol = h.getSolution();
  R.x = sol.col_value;
  R.ineq_marg.assign(sol.row_dual.begin(), sol.row_dual.begin() + nub);
  
  const double tol = std::sqrt(1e-9) * 10;
  bool feas = true;
  for (int j = 0; j < ncol && feas; ++j) {
    const double v = R.x[j];
    if (std::isnan(v) || !(v >= 0.0 - tol)) feas = false;
  }
  for (int i = 0; i < nub && feas; ++i) {
    const double slack = lp.row_upper_[i] - sol.row_value[i];
    if (std::isnan(slack) || slack < -tol) feas = false;
  }
  if (feas) {
    const double con = beq - sol.row_value[nub];
    if (std::isnan(con) || std::fabs(con) > tol) feas = false;
  }
  if (std::isnan(info.objective_function_value)) feas = false;
  R.status = feas ? 0 : 4;
  return R;
}

struct Grid { int size; std::vector<double> A; };  
inline const Grid& grid(int size) {
  static thread_local std::map<int, Grid> cache;
  auto it = cache.find(size);
  if (it != cache.end()) return it->second;
  Grid g; g.size = size; g.A.resize((size_t)D * size);
  const double twopi = 2 * kPi;
  for (int k = 0; k < size; ++k) {
    const double z = twopi * (double)k / (double)size;
    for (int h = 1; h <= K; ++h) {
      const double arg = (double)h * z;
      g.A[(size_t)(h - 1) * size + k] = std::cos(arg);
      g.A[(size_t)(K + h - 1) * size + k] = std::sin(arg);
    }
  }
  return cache.emplace(size, std::move(g)).first->second;
}



inline double lower_polynomial(const double* a) {
  constexpr int N = 131072;
  static const std::vector<double>* tabs = [] {
    auto* t = new std::vector<double>(2 * (size_t)N);
    for (int k = 0; k < N; ++k) { const double th = 2 * kPi * (double)k / N; (*t)[k] = std::cos(th); (*t)[N + k] = std::sin(th); }
    return t;
  }();
  const double* ct = tabs->data(); const double* st = ct + N;
  double vmin = INFINITY;
#pragma omp parallel for reduction(min : vmin) schedule(static)
  for (int j = 0; j < N; ++j) {
    double v = 0;
    for (int h = 1; h <= K; ++h) {
      const int idx = (int)(((long long)h * j) & (N - 1));
      v += a[h - 1] * ct[idx] + a[K + h - 1] * st[idx];
    }
    vmin = std::min(vmin, v);
  }
  double second = 0, asum = 0;
  for (int h = 1; h <= K; ++h) second += (double)h * h * std::hypot(a[h - 1], a[K + h - 1]);
  for (int i = 0; i < D; ++i) asum += std::fabs(a[i]);
  const double step = 2 * kPi / N;
  const double error = second * (step * step) / 8;
  const double rounding = 1e-11 * std::max(1., asum);
  return vmin - error - rounding;
}

inline double disk_support(const double* a, const double* lo, const double* hi, double m) {
  double r[D], tot = 0;
  for (int i = 0; i < D; ++i) r[i] = std::max(a[i] * lo[i], a[i] * hi[i]);
  for (int h = 0; h < K; ++h) tot += std::min(r[h] + r[K + h], m * std::hypot(a[h], a[K + h]));
  return tot;
}

enum class Cone { DiskExcluded, Witness, JointExcluded, Unresolved };

struct GridSolver {
  std::vector<std::pair<std::vector<double>, double>> bank;  
  long long last_iters = 0;                                  
  long long max_lp_iters_resolved = 0;                       

  Cone solve(const Summ& S, double B, double sigma, long long budget) {
    last_iters = 0;
    double lo[D], hi[D], margin;
    const int ex = moment_box(K, S.c, S.s, S.n, B, sigma, lo, hi, &margin);
    const double m = 1 - B;
    if (ex) return Cone::DiskExcluded;
    if (m == 0) {
      const double r = std::sqrt(2 * std::log(51200.0) / S.n);
      double mx = 0;
      for (int h = 0; h < K; ++h) {
        const double v = sigma * ((double)(h + 1) / 4), d = std::exp(-.5 * (v * v));
        mx = std::max({mx, std::fabs(S.c[h] - d), std::fabs(S.s[h])});
      }
      return mx <= r + 1e-12 ? Cone::Witness : Cone::Unresolved;
    }
    for (auto& bk : bank) {
      const double gap = m * bk.second - disk_support(bk.first.data(), lo, hi, m);
      if (gap > 1e-8) return Cone::JointExcluded;
    }
    for (int size : {256, 512, 1024}) {
      const Grid& G = grid(size);
      bool active[2 * D];
      for (int i = 0; i < D; ++i) { active[i] = hi[i] < m; active[D + i] = lo[i] > -m; }
      const int ncol = size + 1;
      std::vector<double> Aub, bub;
      int nub = 0;
      for (int r = 0; r < 2 * D; ++r) {
        if (!active[r]) continue;
        const int row = r % D; const double sg = r < D ? 1.0 : -1.0;
        for (int k = 0; k < size; ++k) Aub.push_back(sg * G.A[(size_t)row * size + k]);
        Aub.push_back(-1.0);
        bub.push_back(r < D ? hi[row] / m : -lo[row] / m);
        ++nub;
      }
      std::vector<double> cost(ncol, 0.0), Aeq(ncol, 1.0);
      cost[size] = 1; Aeq[size] = 0;
      const LPRes res = linprog_highs(ncol, cost, nub, Aub, bub, Aeq, 1.0, budget - last_iters);
      last_iters += res.iters;
      if (res.status != 0) return Cone::Unresolved;  
      max_lp_iters_resolved = std::max(max_lp_iters_resolved, res.iters);
      std::vector<double> w(size);
      double total = 0;
      bool finite = true;
      for (int k = 0; k < size; ++k) { w[k] = std::max(res.x[k], 0.0); finite &= std::isfinite(w[k]); total += w[k]; }
      if (!finite || !std::isfinite(total) || total <= 0) return Cone::Unresolved;
      for (int k = 0; k < size; ++k) w[k] *= m / total;
      double q[D];
      bool qfin = true;
      for (int i = 0; i < D; ++i) {
        double acc = 0;
        for (int k = 0; k < size; ++k) acc += G.A[(size_t)i * size + k] * w[k];
        q[i] = acc; qfin &= std::isfinite(acc);
      }
      if (!qfin) return Cone::Unresolved;
      double viol = 0;
      for (int i = 0; i < D; ++i) viol = std::max({viol, lo[i] - q[i], q[i] - hi[i]});
      const double r = std::sqrt(2 * std::log(51200.0) / S.n);
      double obs = -INFINITY;
      for (int h = 0; h < K; ++h) {
        const double v = sigma * ((double)(h + 1) / 4), d = std::exp(-.5 * (v * v));
        obs = std::max({obs, std::fabs(d * (B + q[h]) - S.c[h]), std::fabs(d * q[K + h] - S.s[h])});
      }
      obs -= r;
      if (viol <= 1e-9 && obs <= 1e-9) return Cone::Witness;
      std::vector<double> lam(2 * D, 0.0);
      for (int r2 = 0, q2 = 0; r2 < 2 * D; ++r2) if (active[r2]) { lam[r2] = std::max(0.0, -res.ineq_marg[q2]); ++q2; }
      std::vector<double> a(D);
      for (int i = 0; i < D; ++i) a[i] = lam[i] - lam[D + i];
      const double lower = lower_polynomial(a.data());
      const double gap = m * lower - disk_support(a.data(), lo, hi, m);
      if (gap > 1e-8) {
        bank.emplace_back(a, lower);
        if (bank.size() > 64) bank.erase(bank.begin(), bank.end() - 64);
        return Cone::JointExcluded;
      }
      if (last_iters >= budget) break;  
    }
    return Cone::Unresolved;
  }
};


struct Screen {
  Summ S;
  GridSolver engine;
  long long iter_first, iter_retry;
  std::map<std::pair<double, double>, std::string> cache;
  std::map<std::pair<double,double>,std::array<int,5>> outcomes;
  
  long long lp_iters = 0, max_call_first = 0, max_call_retry = 0;  
  int n_checks = 0, n_retry = 0, n_unresolved = 0;

  Screen(const Summ& s, long long it1, long long it2) : S(s), iter_first(it1), iter_retry(it2) {}

  
  std::string check(double B, double sigma) {
    const auto key = std::make_pair(B, sigma);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;
    ++n_checks;
    std::array<int,5> state={-1,-1,-1,-1,-1};
    std::string fail;
    double margin; int which;
    const int hit = cdf_screen(S.x, S.f, S.n, B, sigma, &margin, &which);
    state[0]=hit?0:1;
    if (hit) fail = "cdf_inconsistent";
    else {
      double lo[D], hi[D], dm;
      const int ex = moment_box(K, S.c, S.s, S.n, B, sigma, lo, hi, &dm);
      state[1]=ex?0:1;
      if (ex) fail = "moments_inconsistent";
      else {
        const std::string cs = compat_screen(S, B, sigma);
        if (cs == "coordinate_exclusion") {fail = "moments_inconsistent";state[1]=0;}
        else if (cs == "joint_exclusion") {fail = "no_valid_distribution";state[2]=0;state[3]=1;}
        else {
          const Cone r = engine.solve(S, B, sigma, iter_first);
          lp_iters += engine.last_iters;
          if (r != Cone::Unresolved) max_call_first = std::max(max_call_first, engine.last_iters);
          double pm; int pidx;
          const int ph = probit_screen(S.x, S.f, S.n, B, sigma, &pm, &pidx);
          state[4]=ph?0:1;
          state[2]=r==Cone::Witness?1:(r==Cone::Unresolved?-1:0);
          state[3]=r==Cone::Unresolved?0:1;
          if (ph < 0) throw std::runtime_error("probit input error");
          if (r == Cone::DiskExcluded) fail = "moments_inconsistent";
          else if (r == Cone::JointExcluded) fail = "no_valid_distribution";
          else if (ph) fail = "residual_cdf_inconsistent";
          else if (r == Cone::Witness) fail = "";
          else {
            ++n_retry;
            const Cone r2 = engine.solve(S, B, sigma, iter_retry);
            lp_iters += engine.last_iters;
            if (r2 != Cone::Unresolved) max_call_retry = std::max(max_call_retry, engine.last_iters);
            state[2]=r2==Cone::Witness?1:(r2==Cone::Unresolved?-1:0);
            state[3]=r2==Cone::Unresolved?0:1;
            if (r2 == Cone::Witness) fail = "";
            else if (r2 == Cone::JointExcluded || r2 == Cone::DiskExcluded) fail = "no_valid_distribution";
            else { fail = "distribution_check_inconclusive"; ++n_unresolved; }
          }
        }
      }
    }
    outcomes[key]=state;
    cache[key] = fail;
    return fail;
  }
};







inline int check_all(Screen& scr, double B, double sigma) {
  const Summ& S = scr.S;
  int f = 0;
  double margin; int which;
  if (cdf_screen(S.x, S.f, S.n, B, sigma, &margin, &which)) f |= 1;
  double lo[D], hi[D], dm;
  bool moments_inconsistent = moment_box(K, S.c, S.s, S.n, B, sigma, lo, hi, &dm) != 0;
  std::string cs;
  if (!moments_inconsistent) { cs = compat_screen(S, B, sigma); moments_inconsistent = cs == "coordinate_exclusion"; }
  if (moments_inconsistent) f |= 2;
  else if (cs == "joint_exclusion") f |= 4;
  else {
    Cone r = scr.engine.solve(S, B, sigma, scr.iter_first);
    scr.lp_iters += scr.engine.last_iters;
    if (r == Cone::Unresolved) { ++scr.n_retry; r = scr.engine.solve(S, B, sigma, scr.iter_retry); scr.lp_iters += scr.engine.last_iters; }
    if (r == Cone::JointExcluded || r == Cone::DiskExcluded) f |= 4;
    else if (r == Cone::Unresolved) { f |= 8; ++scr.n_unresolved; }
  }
  double pm; int pidx;
  const int ph = probit_screen(S.x, S.f, S.n, B, sigma, &pm, &pidx);
  if (ph < 0) throw std::runtime_error("probit input error");
  if (ph) f |= 16;
  ++scr.n_checks;
  return f;
}


inline bool summ_valid(const Summ& S) {
  for (int h = 0; h < K; ++h) if (!std::isfinite(S.c[h]) || !std::isfinite(S.s[h])) return false;
  for (int j = 0; j < J; ++j) if (!std::isfinite(S.x[j]) || !std::isfinite(S.f[j])) return false;
  for (int j = 1; j < J; ++j) if (S.x[j] < S.x[j - 1] || S.f[j] < S.f[j - 1]) return false;
  return S.n > 0;
}

}  
