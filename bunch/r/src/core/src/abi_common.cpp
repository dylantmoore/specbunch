
#include <cstddef>
#include <string>
#include "sbc_internal.hpp"
#include "specbunch.h"
#ifdef SB_HAVE_BUILD_CONFIG
#include "sb_build_config.h"
#endif

#ifndef SB_VERSION_STRING
#define SB_VERSION_STRING "0.1.0.dev1"
#endif
#ifndef SB_BUILD_FLAGS
#define SB_BUILD_FLAGS "unknown"
#endif

SBC_DEFINE_RESULT_API(sb, SB_API)

extern "C" {
SB_API const char* sb_version(void) { return SB_VERSION_STRING; }
SB_API uint32_t sb_abi_version(void) { return SB_ABI_VERSION; }
SB_API const char* sb_build_info(void) {
  static const std::string s = std::string("compiler=") +
#if defined(__clang__)
      "clang " __clang_version__
#elif defined(__GNUC__)
      "gcc " __VERSION__
#else
      "other"
#endif
      + "; flags=" SB_BUILD_FLAGS + "; openmp=" +
#if defined(_OPENMP)
      "on (" + std::to_string(_OPENMP) + ")"
#else
      "off (single thread)"
#endif
      ;
  return s.c_str();
}
}
