// Standalone detector/guard test; compile outside the source tree with common and omp includes.
#include "sbc_internal.hpp"
#include <cassert>
#include <iostream>
static void env(const char* value) {
#ifdef _WIN32
  _putenv_s("OMP_NUM_THREADS", value ? value : "");
#else
  if(value) setenv("OMP_NUM_THREADS",value,1); else unsetenv("OMP_NUM_THREADS");
#endif
}
int main() {
  env(nullptr);
  const int detected=sbc::default_threads();
  assert(detected>=1);
  { sbc::ThreadGuard guard(0); 
#ifdef _OPENMP
    assert(guard.count()==std::min(detected,omp_get_thread_limit()));
#else
    assert(guard.count()==1);
#endif
#ifdef _OPENMP
    int observed=0;
#pragma omp parallel
    {
#pragma omp single
      observed=omp_get_num_threads();
    }
    assert(observed==guard.count());
#endif
  }
#ifdef __APPLE__
  int performance=0, physical=0; size_t size=sizeof(int);
  if(sysctlbyname("hw.perflevel0.physicalcpu",&performance,&size,nullptr,0)==0) {
    assert(detected==performance);
    size=sizeof(int); sysctlbyname("hw.physicalcpu",&physical,&size,nullptr,0);
    std::cout << "Mac performance=" << performance << " physical=" << physical << '\n';
  }
#endif
  for(const char* bad : {"", "0", "-1", "garbage", "999999999999999999999"}) {
    env(bad); assert(sbc::default_threads()==detected);
  }
  env("3"); assert(sbc::default_threads()==3);
  env("3,2"); assert(sbc::default_threads()==3);
  for(int requested : {0,-1,1,3,16}) {
    int before=omp_get_max_threads();
    { sbc::ThreadGuard guard(requested);
      int expected=requested>0 ? requested : 3;
#ifdef _OPENMP
      assert(guard.count()==std::min(expected,omp_get_thread_limit()));
      assert(omp_get_max_threads()==expected);
      int observed=0;
#pragma omp parallel
      {
#pragma omp single
        observed=omp_get_num_threads();
      }
      assert(observed==guard.count());
#else
      (void)expected;
      assert(guard.count()==1);
#endif
    }
    assert(omp_get_max_threads()==before);
  }
  env(nullptr);
  std::cout << "detector, override, explicit 1/3/16, restoration PASS\n";
}
