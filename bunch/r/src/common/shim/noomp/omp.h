/* omp.h stand-in for builds without OpenMP (sbcommon). The build adds this directory to the include path ONLY when no
 * OpenMP runtime is found; the pragmas are then ignored (the build uses -fopenmp-simd, so `omp simd` reductions keep
 * the same lane structure as an OpenMP build) and every parallel region runs on one thread. The cores use fixed-order
 * reductions, so results equal those of an OpenMP build at any thread count. */
#ifndef SBC_NOOMP_OMP_H
#define SBC_NOOMP_OMP_H
#ifdef __cplusplus
extern "C" {
#endif
static inline int omp_get_max_threads(void) { return 1; }
static inline int omp_get_num_threads(void) { return 1; }
static inline int omp_get_thread_num(void) { return 0; }
static inline int omp_get_num_procs(void) { return 1; }
static inline void omp_set_num_threads(int n) { (void)n; }
static inline void omp_set_max_active_levels(int n) { (void)n; }
static inline int omp_in_parallel(void) { return 0; }
#ifdef __cplusplus
}
#endif
#endif
