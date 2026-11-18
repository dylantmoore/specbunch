



#pragma once
#include <algorithm>
#include <cmath>
#include <vector>
#include "distribution_screen.hpp"

namespace sbdistribution {


inline void summarize(const std::vector<double>& y, double kink, double* rec) {
  constexpr int K = 128, J = 66;
  const size_t n = y.size();
  long double C[K] = {0}, S[K] = {0};
  constexpr size_t BLK = 1024;
  double c[K], s[K];
  for (size_t i0 = 0; i0 < n; i0 += BLK) {
    size_t i1 = std::min(n, i0 + BLK);
    std::fill(c, c + K, 0.0); std::fill(s, s + K, 0.0);
    for (size_t i = i0; i < i1; ++i) {
      const double z = y[i] - kink;
      const double rc = std::cos(0.25 * z), rs = std::sin(0.25 * z);
      for (int h0 = 1; h0 <= K; h0 += 8) {
        double tc = std::cos(0.25 * h0 * z), ts = std::sin(0.25 * h0 * z);
        for (int h = h0; h < h0 + 8 && h <= K; ++h) {
          c[h - 1] += tc; s[h - 1] += ts;
          double nc = tc * rc - ts * rs, ns = ts * rc + tc * rs; tc = nc; ts = ns;
        }
      }
    }
    for (int h = 0; h < K; ++h) { C[h] += c[h]; S[h] += s[h]; }
  }
  rec[0] = (double)n;
  for (int h = 0; h < K; ++h) { rec[1 + h] = (double)(C[h] / n); rec[1 + K + h] = (double)(S[h] / n); }
  std::vector<double> v(y); std::sort(v.begin(), v.end());
  std::vector<double> x;
  for (int j = 0; j <= 64; ++j) x.push_back(v[(long long)j * (n - 1) / 64]);
  x.push_back(kink); std::sort(x.begin(), x.end());
  for (int j = 0; j < J; ++j) {
    rec[1 + 2 * K + j] = x[j] - kink;
    rec[1 + 2 * K + J + j] = double(std::upper_bound(v.begin(), v.end(), x[j]) - v.begin()) / n;
  }
}

inline void summarize_w(const std::vector<double>& y, const std::vector<double>& w, double* rec) {
  constexpr int K = 128, J = 66;
  const size_t n = y.size();
  long double C[K] = {0}, S[K] = {0}, Wt = 0, W2 = 0;
  constexpr size_t BLK = 1024;
  double c[K], s[K];
  for (size_t i0 = 0; i0 < n; i0 += BLK) {
    size_t i1 = std::min(n, i0 + BLK);
    std::fill(c, c + K, 0.0); std::fill(s, s + K, 0.0);
    double wb = 0, w2b = 0;
    for (size_t i = i0; i < i1; ++i) {
      const double z = y[i], wi = w[i];
      wb += wi; w2b += wi * wi;
      const double rc = std::cos(0.25 * z), rs = std::sin(0.25 * z);
      for (int h0 = 1; h0 <= K; h0 += 8) {
        double tc = std::cos(0.25 * h0 * z), ts = std::sin(0.25 * h0 * z);
        for (int h = h0; h < h0 + 8 && h <= K; ++h) {
          c[h - 1] += wi * tc; s[h - 1] += wi * ts;
          double nc = tc * rc - ts * rs, ns = ts * rc + tc * rs; tc = nc; ts = ns;
        }
      }
    }
    for (int h = 0; h < K; ++h) { C[h] += c[h]; S[h] += s[h]; }
    Wt += wb; W2 += w2b;
  }
  rec[0] = std::floor((double)(Wt * Wt / W2));
  for (int h = 0; h < K; ++h) { rec[1 + h] = (double)(C[h] / Wt); rec[1 + K + h] = (double)(S[h] / Wt); }
  std::vector<size_t> ord(n);
  for (size_t i = 0; i < n; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return y[a] < y[b] || (y[a] == y[b] && a < b); });
  std::vector<double> v(n), cw(n);
  long double acc = 0;
  for (size_t i = 0; i < n; ++i) { v[i] = y[ord[i]]; acc += w[ord[i]]; cw[i] = (double)(acc / Wt); }
  std::vector<double> x;
  x.push_back(v[0]);
  for (int j = 1; j < 64; ++j) {
    const double q = j / 64.0;
    size_t k = (size_t)(std::lower_bound(cw.begin(), cw.end(), q) - cw.begin());
    x.push_back(v[std::min(k, n - 1)]);
  }
  x.push_back(v[n - 1]);
  x.push_back(0.0); std::sort(x.begin(), x.end());
  for (int j = 0; j < J; ++j) {
    rec[1 + 2 * K + j] = x[j];
    size_t k = (size_t)(std::upper_bound(v.begin(), v.end(), x[j]) - v.begin());
    rec[1 + 2 * K + J + j] = k == 0 ? 0.0 : std::min(1.0, cw[k - 1]);
  }
}


inline Summ to_summ(const double* rec) {
  Summ S;
  S.n = (int)rec[0];
  for (int h = 0; h < 128; ++h) { S.c[h] = rec[1 + h]; S.s[h] = rec[129 + h]; }
  for (int j = 0; j < 66; ++j) { S.x[j] = rec[257 + j]; S.f[j] = rec[323 + j]; }
  return S;
}

}  
