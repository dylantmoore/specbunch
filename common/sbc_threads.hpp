// Shared default thread policy. No estimator or wrapper owns CPU detection.
#pragma once
#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#if defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#elif defined(__linux__)
#include <dirent.h>
#include <sched.h>
#elif defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace sbc {
inline int positive_env_threads() {
  const char* value = std::getenv("OMP_NUM_THREADS");
  if (!value) return 0;
  char* end = nullptr;
  errno = 0;
  long n = std::strtol(value, &end, 10);
  if (end == value || errno || n < 1 || n > INT_MAX) return 0;
  while (*end == ' ' || *end == '\t') ++end;
  // OpenMP permits a comma-separated list for nested regions: use the first level.
  return (*end == '\0' || *end == ',') ? static_cast<int>(n) : 0;
}
#if defined(__linux__)
inline std::set<int> cpu_list(const std::string& s) {
  std::set<int> out;
  std::istringstream in(s); std::string part;
  while (std::getline(in, part, ',')) {
    std::istringstream p(part); int a, b;
    if (!(p >> a) || a < 0) continue;
    b = a;
    if (p.peek() == '-') { p.get(); if (!(p >> b)) continue; }
    if (b < a || b > 1048576) continue;
    for (int i=a; i<=b; ++i) out.insert(i);
  }
  return out;
}
inline std::string read_line(const std::string& path) {
  std::ifstream f(path); std::string s; std::getline(f,s); return s;
}
inline int linux_threads(int fallback) {
  std::set<int> allowed;
  // Grow the affinity buffer for machines with more than CPU_SETSIZE CPUs.
  for (int n=128; n<=1048576; n*=2) {
    cpu_set_t* mask=CPU_ALLOC(n);
    if (!mask) break;
    size_t bytes=CPU_ALLOC_SIZE(n); CPU_ZERO_S(bytes,mask);
    int rc=sched_getaffinity(0,bytes,mask), error=errno;
    if (!rc) for(int i=0;i<n;++i) if(CPU_ISSET_S(i,bytes,mask)) allowed.insert(i);
    CPU_FREE(mask);
    if (!rc || error != EINVAL) break;
  }
  int cap=allowed.empty() ? fallback : static_cast<int>(allowed.size());
  // Count machine cores, then cap by the calling process CPU allowance.
  auto performance=cpu_list(read_line("/sys/devices/cpu_core/cpus"));
  auto online=cpu_list(read_line("/sys/devices/system/cpu/online"));
  std::set<std::pair<int,int>> all, fast;
  if (DIR* dir=opendir("/sys/devices/system/cpu")) {
    while (auto* ent=readdir(dir)) {
      std::string name=ent->d_name;
      if(name.compare(0,3,"cpu") || name.size()==3) continue;
      char* end=nullptr; long cpu=std::strtol(name.c_str()+3,&end,10);
      if(*end || cpu<0 || cpu>INT_MAX) continue;
      if(!online.empty() && !online.count(static_cast<int>(cpu))) continue;
      std::string base="/sys/devices/system/cpu/"+name+"/topology/";
      int package=-1, core=-1;
      std::istringstream(read_line(base+"physical_package_id")) >> package;
      std::istringstream(read_line(base+"core_id")) >> core;
      if(package<0 || core<0) continue;
      all.emplace(package,core);
      if(performance.count(static_cast<int>(cpu))) fast.emplace(package,core);
    }
    closedir(dir);
  }
  // cgroup v2 quota, including ancestor quotas. Missing files are harmless.
  std::ifstream groups("/proc/self/cgroup"); std::string line;
  while(std::getline(groups,line)) if(line.rfind("0::/",0)==0) {
    std::string relative=line.substr(3);
    if(relative.find("..") != std::string::npos) break;
    std::string path="/sys/fs/cgroup"+(relative=="/" ? "" : relative);
    for (;;) {
      std::istringstream q(read_line(path+"/cpu.max"));
      long long quota=0, period=0;
      if((q>>quota>>period) && quota>0 && period>0)
        cap=std::min(cap,static_cast<int>(std::min<long long>(INT_MAX,std::max(1LL,quota/period))));
      if(path=="/sys/fs/cgroup") break;
      auto slash=path.rfind('/');
      if(slash<14 || slash==std::string::npos) break;
      path.resize(slash);
    }
    break;
  }
  int cores=static_cast<int>(fast.empty() ? all.size() : fast.size());
  return std::max(1,std::min(cores ? cores : fallback,cap));
}
#endif
inline int detected_threads() {
  int fallback=static_cast<int>(std::thread::hardware_concurrency());
  fallback=std::max(1,fallback);
#if defined(__APPLE__)
  auto query=[](const char* key) { int n=0; size_t size=sizeof(n);
    return sysctlbyname(key,&n,&size,nullptr,0)==0 && n>0 ? n : 0; };
  int cores=query("hw.perflevel0.physicalcpu");
  if(!cores) cores=query("hw.physicalcpu");
  int available=query("hw.activecpu");
  return std::max(1,std::min(cores ? cores : fallback,available ? available : fallback));
#elif defined(__linux__)
  return linux_threads(fallback);
#elif defined(_WIN32)
  DWORD_PTR process=0, system=0;
  int cap=fallback;
  bool affinity=GetProcessAffinityMask(GetCurrentProcess(),&process,&system)!=0 && process;
  if(affinity) { cap=0; for(auto m=process;m;m>>=1) cap+=static_cast<int>(m&1); }
  else {
    // Multi-group processes may have a zero legacy affinity mask.
    USHORT count=GetActiveProcessorGroupCount();
    std::vector<USHORT> groups(count);
    if(count && GetProcessGroupAffinity(GetCurrentProcess(),&count,groups.data())) {
      int total=0;
      for(USHORT i=0;i<count;++i) total+=static_cast<int>(GetActiveProcessorCount(groups[i]));
      if(total>0) cap=total;
    }
  }
  DWORD bytes=0;
  GetLogicalProcessorInformationEx(RelationProcessorCore,nullptr,&bytes);
  std::vector<unsigned char> storage(bytes);
  int best=-1, cores=0;
  if(bytes && GetLogicalProcessorInformationEx(RelationProcessorCore,
      reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(storage.data()),&bytes)) {
    for(size_t offset=0;offset<bytes;) {
      auto* info=reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(storage.data()+offset);
      if(!info->Size || offset+info->Size>bytes) break;
      if(info->Relationship==RelationProcessorCore) {
        int efficiency=info->Processor.EfficiencyClass;
        if(efficiency>best) { best=efficiency; cores=0; }
        if(efficiency==best) ++cores;
      }
      offset+=info->Size;
    }
  }
  return std::max(1,std::min(cores ? cores : fallback,cap));
#else
  return fallback;
#endif
}
// Explicit counts and user OMP_NUM_THREADS are deliberate overrides, uncapped.
inline int default_threads() noexcept {
  if(int n=positive_env_threads()) return n;
  try { return detected_threads(); } catch (...) { return 1; }
}
} // namespace sbc
