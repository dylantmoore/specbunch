// Test-only bridge linked to the unmodified Homebrew dynamic runtime.
#include <omp.h>
extern "C" __attribute__((visibility("default"))) int foreign_parallel(int n) {
  int observed = 0;
#pragma omp parallel num_threads(n) reduction(+:observed)
  observed += 1;
  return observed;
}
