// sbc_internal.hpp: C++ side of the sbcommon ABI conventions (sb_common.h). Compiled into each package as source.
// Nothing here crosses the C boundary except through the prefixed extern "C" functions that SBC_DEFINE_RESULT_API
// generates in one translation unit of each package.
#pragma once
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <map>
#include <new>
#include <string>
#include <vector>
#include <omp.h>  // the real one, or common/shim/noomp/omp.h in a build without OpenMP
#include "sb_common.h"
#include "sbc_threads.hpp"

namespace sbc {

// ---------------------------------------------------------------- result handle
struct Result {
  int status = SBC_STATUS_ESTIMATE;
  std::map<std::string, double> scalars;
  std::map<std::string, std::vector<double>> vectors;
  std::map<std::string, std::string> strings;
  void set(const std::string& k, double v) { scalars[k] = v; }
  void setv(const std::string& k, std::vector<double> v) { vectors[k] = std::move(v); }
  void sets(const std::string& k, const std::string& v) { strings[k] = v; }
  // newline-separated list of names, for discovery from any language
  std::string names_scalar, names_vector, names_string;
  void finish_names() {
    names_scalar.clear(); names_vector.clear(); names_string.clear();
    for (auto& kv : scalars) names_scalar += (names_scalar.empty() ? "" : "\n") + kv.first;
    for (auto& kv : vectors) names_vector += (names_vector.empty() ? "" : "\n") + kv.first;
    for (auto& kv : strings) names_string += (names_string.empty() ? "" : "\n") + kv.first;
  }
};

// ---------------------------------------------------------------- thread-local last error
inline std::string& last_error_ref() { static thread_local std::string e; return e; }
inline void clear_error() { last_error_ref().clear(); }
inline int fail(int code, const std::string& msg) { last_error_ref() = msg; return code; }

// ---------------------------------------------------------------- OpenMP thread count: set for the call, restore after
class ThreadGuard {
 public:
  explicit ThreadGuard(int threads) : prev_(omp_get_max_threads()), count_(1) {
#ifdef _OPENMP
    omp_set_num_threads(threads > 0 ? threads : default_threads());
    // Observe a real team, including runtime thread limits/nested serialization.
    // Do not mistake the requested maximum for the number of workers available.
#pragma omp parallel
    {
#pragma omp single
      count_ = omp_get_num_threads();
    }
#else
    (void)threads;  // noomp / -fopenmp-simd: parallel regions execute serially
#endif
  }
  int count() const { return count_; }
  ~ThreadGuard() { omp_set_num_threads(prev_); }
  ThreadGuard(const ThreadGuard&) = delete;
  ThreadGuard& operator=(const ThreadGuard&) = delete;
 private:
  int prev_, count_;
};

// ---------------------------------------------------------------- options struct header check
// Every options struct begins with {uint32 struct_size, uint32 abi_version}. min_size = offset of the end of the last
// field the core requires; fields beyond struct_size take their defaults.
inline int check_opts_header(const void* opts, uint32_t min_size, uint32_t max_version) {
  if (!opts) return fail(SBC_E_NULL, "options pointer is NULL");
  uint32_t sz, ver;
  std::memcpy(&sz, opts, 4); std::memcpy(&ver, (const char*)opts + 4, 4);
  if (sz < min_size) return fail(SBC_E_OPTS, "options struct_size " + std::to_string(sz) + " < required " + std::to_string(min_size));
  if (ver == 0 || ver > max_version)
    return fail(SBC_E_OPTS, "options abi_version " + std::to_string(ver) + " not supported (max " + std::to_string(max_version) + ")");
  return SBC_OK;
}

// Run f() with every exception mapped to an error code (no exception crosses the C boundary).
template <class F>
inline int guarded(F&& f) {
  clear_error();
  try { return f(); }
  catch (const std::bad_alloc&) { return fail(SBC_E_INTERNAL, "out of memory"); }
  catch (const std::exception& e) { return fail(SBC_E_INTERNAL, std::string("internal error: ") + e.what()); }
  catch (...) { return fail(SBC_E_INTERNAL, "internal error (unknown exception)"); }
}

}  // namespace sbc

// Defines the prefixed result getters in one translation unit of a package, e.g. SBC_DEFINE_RESULT_API(sb, SB_API).
// The result type is the package's opaque struct PREFIX##_result, which is an sbc::Result.
#define SBC_DEFINE_RESULT_API(P, API)                                                                                   \
  extern "C" {                                                                                                          \
  API int P##_default_threads(void) { return sbc::default_threads(); } \
  API const char* P##_last_error(void) { return sbc::last_error_ref().c_str(); }                                        \
  API int P##_result_status(const P##_result* r) { return r ? reinterpret_cast<const sbc::Result*>(r)->status : -1; }  \
  API int P##_result_scalar(const P##_result* r, const char* name, double* out) {                                       \
    if (!r || !name || !out) return sbc::fail(SBC_E_NULL, "NULL argument");                                             \
    auto& m = reinterpret_cast<const sbc::Result*>(r)->scalars;                                                         \
    auto it = m.find(name);                                                                                             \
    if (it == m.end()) { *out = std::nan(""); return sbc::fail(SBC_E_NOTFOUND, std::string("no scalar ") + name); }     \
    *out = it->second; return SBC_OK;                                                                                   \
  }                                                                                                                     \
  API const double* P##_result_vector(const P##_result* r, const char* name, size_t* len) {                             \
    if (len) *len = 0;                                                                                                  \
    if (!r || !name) { sbc::fail(SBC_E_NULL, "NULL argument"); return nullptr; }                                        \
    auto& m = reinterpret_cast<const sbc::Result*>(r)->vectors;                                                         \
    auto it = m.find(name);                                                                                             \
    if (it == m.end()) { sbc::fail(SBC_E_NOTFOUND, std::string("no vector ") + name); return nullptr; }                 \
    if (len) *len = it->second.size();                                                                                  \
    return it->second.data();                                                                                           \
  }                                                                                                                     \
  API const char* P##_result_string(const P##_result* r, const char* name) {                                            \
    if (!r || !name) { sbc::fail(SBC_E_NULL, "NULL argument"); return nullptr; }                                        \
    auto& m = reinterpret_cast<const sbc::Result*>(r)->strings;                                                         \
    auto it = m.find(name);                                                                                             \
    if (it == m.end()) { sbc::fail(SBC_E_NOTFOUND, std::string("no string ") + name); return nullptr; }                 \
    return it->second.c_str();                                                                                          \
  }                                                                                                                     \
  API const char* P##_result_names(const P##_result* r, int kind) {                                                     \
    if (!r) return "";                                                                                                  \
    auto* q = reinterpret_cast<const sbc::Result*>(r);                                                                  \
    return kind == 0 ? q->names_scalar.c_str() : kind == 1 ? q->names_vector.c_str() : q->names_string.c_str();         \
  }                                                                                                                     \
  API void P##_result_free(P##_result* r) { delete reinterpret_cast<sbc::Result*>(r); }                                 \
  }
