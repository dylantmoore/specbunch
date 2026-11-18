// Portable timing compatibility for unchanged vendor sources. Timing never selects an estimate.
#pragma once
#include <ctime>
struct sb_rusage_timeval { long tv_sec; long tv_usec; };
struct rusage { sb_rusage_timeval ru_utime, ru_stime; };
#define RUSAGE_SELF 0
static inline int getrusage(int, rusage* r) {
  const double c = double(std::clock()) / CLOCKS_PER_SEC;
  r->ru_utime.tv_sec = long(c);
  r->ru_utime.tv_usec = long((c - double(long(c))) * 1e6);
  r->ru_stime = {0, 0};
  return 0;
}
