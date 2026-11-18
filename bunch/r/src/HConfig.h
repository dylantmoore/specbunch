#ifndef HCONFIG_H_
#define HCONFIG_H_
#define FAST_BUILD
#define CMAKE_BUILD_TYPE "Release"
#define HIGHS_GITHASH "n/a"
#define HIGHS_VERSION_MAJOR 1
#define HIGHS_VERSION_MINOR 8
#define HIGHS_VERSION_PATCH 0
#define HIGHS_DIR "."
#if defined(__GNUC__) || defined(__clang__)
#define HIGHS_HAVE_BUILTIN_CLZ
#endif
#endif
