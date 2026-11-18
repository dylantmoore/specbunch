/* Adapter-only changes: reference files remain byte-identical. */
#ifndef SPECTRAL_R_PORTABILITY_H
#define SPECTRAL_R_PORTABILITY_H
#include "r_build_config.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#ifdef __cplusplus
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
namespace spectral_r {
struct nullbuf : std::streambuf { int overflow(int c) override { return c; } };
inline nullbuf sink;
inline std::ostream quiet(&sink);
[[noreturn]] inline void stop(int) { throw std::runtime_error("Disabled standalone program exit"); }
[[noreturn]] inline void abort_call() { throw std::runtime_error("Native solver abort"); }
inline int random_call() { throw std::runtime_error("Unexpected C rand call"); }
inline int quiet_printf(const char*, ...) { return 0; }
inline int quiet_fprintf(FILE*, const char*, ...) { return 0; }
inline int quiet_fputs(const char*, FILE*) { return 0; }
inline int quiet_vfprintf(FILE*, const char*, va_list) { return 0; }
}
namespace std { using spectral_r::quiet; using spectral_r::stop; using spectral_r::abort_call; using spectral_r::random_call; using spectral_r::quiet_printf; using spectral_r::quiet_fprintf; using spectral_r::quiet_vfprintf; using spectral_r::quiet_fputs; }
using spectral_r::stop; using spectral_r::abort_call; using spectral_r::random_call;
using spectral_r::quiet_printf; using spectral_r::quiet_fprintf; using spectral_r::quiet_vfprintf; using spectral_r::quiet_fputs;
#define cout quiet
#define cerr quiet
#define exit stop
#define abort abort_call
#define rand random_call
#define printf quiet_printf
#define fprintf quiet_fprintf
#define vfprintf quiet_vfprintf
#define fputs quiet_fputs
#else
static inline int spectral_quiet_printf(const char* fmt, ...) { (void)fmt; return 0; }
static inline int spectral_quiet_fprintf(FILE* stream, const char* fmt, ...) { (void)stream; (void)fmt; return 0; }
static inline int spectral_quiet_fputs(const char* s, FILE* stream) { (void)s; (void)stream; return 0; }
void spectral_r_stop(int code);
#define exit spectral_r_stop
#define fprintf spectral_quiet_fprintf
#define fputs spectral_quiet_fputs
#define printf spectral_quiet_printf
#endif
#undef stdout
#undef stderr
#define stdout ((FILE*)0)
#define stderr ((FILE*)0)
#endif
