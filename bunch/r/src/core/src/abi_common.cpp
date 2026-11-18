
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

extern "C" SB_API int sb_mass(const sb_data* data,const sb_opts_mass* opts,int32_t method,int32_t friction,sb_result** out) {
  return sbc::guarded([&] {
    if(out)*out=nullptr;
    if(method!=SB_METHOD_COSINE && method!=SB_METHOD_JOINT)return sbc::fail(SBC_E_INPUT,"method must be cosine or joint");
    if(friction!=SB_FRICTION_GAUSSIAN && friction!=SB_FRICTION_VOIGT)return sbc::fail(SBC_E_INPUT,"friction must be gaussian or voigt");
    if(method==SB_METHOD_COSINE)return friction==SB_FRICTION_GAUSSIAN?sb_mass_gaussian(data,opts,out):sb_mass_voigt(data,opts,out);
    return friction==SB_FRICTION_GAUSSIAN?sb_mass_gaussian_jump(data,opts,out):sb_mass_joint_voigt(data,opts,out);
  });
}
