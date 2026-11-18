/* sb_common.h: ABI conventions shared by the spectral-bunching packages (sbcommon).
 *
 * Pure C (C99). This header holds conventions only: no function is declared here. Each package declares its own
 * prefixed entry points (sb_ for specbunch, srd_ for specrd) in its public header, so both libraries can be loaded
 * in one process. sbcommon is compiled into each package's binary as source; it is never a runtime library.
 *
 * Conventions every package follows:
 *  1. Options structs start with two fields, `uint32_t struct_size; uint32_t abi_version;`. Callers set them with the
 *     package's *_init() function (or by hand: struct_size = sizeof(struct), abi_version = the header's version). The
 *     core reads only the fields that fit in struct_size, so fields can be appended without breaking old callers.
 *  2. Results come back through an opaque handle with named getters (scalar, vector, string), a status code
 *     (SBC_STATUS_ESTIMATE = 0, > 0 = refused, with a reason string) and a free function.
 *  3. Calls return an int error code (SBC_OK = 0, else one of SBC_E_*). The message for the calling thread's last error
 *     is available from the package's *_last_error() (thread-local; empty string after a successful call).
 *  4. No C++ exception, std:: type or ownership of caller memory crosses the boundary.
 *  5. Thread count is an argument (`threads`, <= 0 = OMP_NUM_THREADS or detected performance/physical cores). The core sets the OpenMP thread count for the call and
 *     restores the caller's setting afterwards. Results are bitwise identical at any thread count.
 *  6. Missing values: rows with a non-finite value (or, if weights are given, a non-finite or non-positive weight) are
 *     dropped inside the core, so every language drops the same rows. The getters report the n actually used.
 */
#ifndef SB_COMMON_H
#define SB_COMMON_H

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32) || defined(__CYGWIN__)
#  if defined(SBC_BUILDING)
#    define SBC_API __declspec(dllexport)
#  else
#    define SBC_API __declspec(dllimport)
#  endif
#else
#  define SBC_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
#  define SBC_EXTERN_C_BEGIN extern "C" {
#  define SBC_EXTERN_C_END }
#else
#  define SBC_EXTERN_C_BEGIN
#  define SBC_EXTERN_C_END
#endif

/* Error codes returned by every entry point. */
#define SBC_OK 0
#define SBC_E_NULL 1        /* a required pointer is NULL */
#define SBC_E_OPTS 2        /* options struct: struct_size too small or abi_version unsupported */
#define SBC_E_INPUT 3       /* bad input (e.g. n == 0, an unsupported option value) */
#define SBC_E_INTERNAL 4    /* an internal error (caught exception, allocation failure) */
#define SBC_E_NOTFOUND 5    /* a getter's name is not in the result */

/* Result status (the estimator's answer, not a call error). */
#define SBC_STATUS_ESTIMATE 0   /* an estimate is returned */
/* > 0: refused; the package header lists its codes; the reason string names the check that failed */

/* The options-struct header every options struct starts with. */
#define SBC_OPTS_HEADER uint32_t struct_size; uint32_t abi_version;

#endif /* SB_COMMON_H */
