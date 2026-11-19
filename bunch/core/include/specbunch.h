/* Public C interface; named results are documented in docs/REFERENCE.md. */
#ifndef SPECBUNCH_H
#define SPECBUNCH_H
#include "sb_common.h"
#define SB_API SBC_API
#define SB_ABI_VERSION 1u
SBC_EXTERN_C_BEGIN
#define SB_POLICY_TIERED 0
#define SB_POLICY_FIRST 2
#define SB_REFUSED_INSUFFICIENT_DATA 1
#define SB_REFUSED_NO_FREQUENCY_BAND 2
#define SB_REFUSED_NO_COHERENT_WINDOW 3
typedef struct sb_result sb_result;
typedef struct sb_data {size_t n;const double *z;const double *w;} sb_data;
typedef struct sb_data_xy {size_t n;const double *z;const double *y;const double *w;} sb_data_xy;
typedef struct sb_opts_mass {
 SBC_OPTS_HEADER
 int32_t threads;
 int32_t check_all_windows; /* 0: selection walk; 2: all coherent windows */
 uint64_t seed;
 int32_t distribution_checks; /* private reserved switch; initialized to 1 */
 int32_t policy; /* SB_POLICY_TIERED or SB_POLICY_FIRST */
 int64_t lp_iterations;
 int64_t lp_retry_iterations;
} sb_opts_mass;
typedef struct sb_opts_bmean {
 SBC_OPTS_HEADER
 int32_t threads;
 int32_t reserved_band;
 int32_t reserved_visibility;
 int32_t R; /* number of fit-test draws */
 uint64_t seed;
} sb_opts_bmean;
typedef struct sb_opts_npmle {
 SBC_OPTS_HEADER
 int32_t threads, core_intervals, grid_points, exact_threshold;
 int32_t max_iterations, lambda_points, refinement_points, reserved; /* reserved: 0 cubic (default), 1 linear */
 double B, sigma, a, kink, histogram_width, lambda_min, lambda_max;
 double quadrature_tolerance, optimality_tolerance;
} sb_opts_npmle;
/* Postestimation options: additive ABI; prior structures are unchanged. */
typedef struct sb_opts_nonbuncher_mean {
 SBC_OPTS_HEADER
 int32_t threads, degree;
 double B, sigma, a, kink, mu_b, se_mu_b, bandwidth, bin_width;
} sb_opts_nonbuncher_mean;
typedef struct sb_opts_eti {
 SBC_OPTS_HEADER
 int32_t threads, integration_intervals;
 double B, sigma, a, kink, window_lo, window_hi, net_below, net_above;
} sb_opts_eti;
SB_API void sb_opts_nonbuncher_mean_init(sb_opts_nonbuncher_mean*);
SB_API void sb_opts_eti_init(sb_opts_eti*);
SB_API int sb_nonbuncher_mean(const sb_data_xy*,const sb_opts_nonbuncher_mean*,sb_result**);
SB_API int sb_eti(const sb_data*,const sb_opts_eti*,sb_result**);
/* ETI from explicit latent density limits; no sample or frequency integration. */
SB_API int sb_eti_shoulders(const sb_opts_eti*,double,double,sb_result**);
SB_API void sb_opts_npmle_init(sb_opts_npmle*);
SB_API int sb_npmle(const sb_data*, const sb_opts_npmle*, const double*, size_t, sb_result**);
SB_API const char* sb_version(void);
SB_API uint32_t sb_abi_version(void);
SB_API const char* sb_build_info(void);
SB_API int sb_default_threads(void);
SB_API const char* sb_last_error(void);
SB_API void sb_opts_mass_init(sb_opts_mass*);
SB_API void sb_opts_bmean_init(sb_opts_bmean*);
#define SB_METHOD_COSINE 0
#define SB_METHOD_JOINT 1
#define SB_FRICTION_GAUSSIAN 0
#define SB_FRICTION_VOIGT 1
/* Additive ABI: option layouts and prior entry points remain compatible. */
SB_API int sb_mass(const sb_data*,const sb_opts_mass*,int32_t method,int32_t friction,sb_result**);
SB_API int sb_mass_joint_voigt(const sb_data*,const sb_opts_mass*,sb_result**);
SB_API int sb_mass_gaussian(const sb_data*,const sb_opts_mass*,sb_result**);
SB_API int sb_mass_voigt(const sb_data*,const sb_opts_mass*,sb_result**);
SB_API int sb_mass_gaussian_jump(const sb_data*,const sb_opts_mass*,sb_result**);
SB_API int sb_bmean(const sb_data_xy*,const sb_opts_bmean*,sb_result**);
SB_API int sb_result_status(const sb_result*);
SB_API int sb_result_scalar(const sb_result*,const char*,double*);
SB_API const double* sb_result_vector(const sb_result*,const char*,size_t*);
SB_API const char* sb_result_string(const sb_result*,const char*);
SB_API const char* sb_result_names(const sb_result*,int);
SB_API void sb_result_free(sb_result*);
SBC_EXTERN_C_END
#endif
