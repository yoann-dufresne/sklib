cmake_minimum_required(VERSION 3.22.0)
project(sklib VERSION 0.13.2 LANGUAGES C CXX)

set(PROJECT_DESCRIPTION "Sorted k-mer library — compact sort and query of k-mers")

# set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

cmake_policy(SET CMP0135 NEW)

option(PORTABLE_BUILD "Generic x86_64 build" OFF)
option(WITH_APP "Compile CLI applications" ON)
option(WITH_TESTS "Compile tests" ON)
option(STATIC_BUILD "Static build" OFF)
option(WITH_UNION_BENCH "Compile the set_union mono-thread benchmark harness" OFF)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE="Release")
endif()

set(CMAKE_CXX_STANDARD 23)

# The widest packed-integer backend (kuint256) is `unsigned _BitInt(256)`, used to index large k-mers
# (2*(2k-m) up to 512 bits). C++ _BitInt(256) needs Clang >= 16 or GCC >= 15 (GCC <= 14 supports
# _BitInt only in C). Probe the actual compiler and fail configuration early with clear guidance
# rather than letting a stale compiler hit a cryptic parse error in the headers.
include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_FLAGS "-std=c++23")
check_cxx_source_compiles(
  "unsigned _BitInt(256) f(unsigned _BitInt(256) x){return x*x>>3;} int main(){return (int)f(5);}"
  SKLIB_HAS_BITINT256)
unset(CMAKE_REQUIRED_FLAGS)
if(NOT SKLIB_HAS_BITINT256)
  message(FATAL_ERROR
    "sklib needs C++ _BitInt(256) support (Clang >= 16 or GCC >= 15) for its 256-bit k-mer backend, "
    "but ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} lacks it (note: g++ <= 14 supports "
    "_BitInt only in C, not C++). On Ubuntu 24.04 install Clang and re-run cmake with it, e.g. "
    "`apt install clang-18` then -DCMAKE_CXX_COMPILER=clang++-18.")
endif()

# Interprocedural optimization (LTO). The lib/app targets carry `if (ipo_supported)`
# INTERPROCEDURAL_OPTIMIZATION_RELEASE blocks, but `ipo_supported` was never defined — so LTO never
# fired. Probe it here, once, after the compiler is known; the targets enable it (Release only) so
# cross-TU calls — notably km::chaining::colinear_chaining (lib/src/algorithms/ColinearChaining.cpp),
# invoked per column in the phase-2 compaction hotspot — can inline. Scoped to our targets (not a
# global CMAKE_INTERPROCEDURAL_OPTIMIZATION) so the third-party deps with an older policy (CMP0069)
# are untouched, and Debug/Profile stay LTO-free/fast. Measured: ~+3% end-to-end construct at -t8
# (chr21/celegans, k31), byte-identical output across all widths (k21/k31/k63).
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error LANGUAGES CXX)
if(ipo_supported)
  message(STATUS "sklib: IPO/LTO available — enabled for Release targets")
else()
  message(STATUS "sklib: IPO/LTO not supported by toolchain: ${ipo_error}")
endif()

add_library(build_type_flags INTERFACE)
target_compile_options(build_type_flags INTERFACE
  $<$<CONFIG:Debug>:-O -g>
  $<$<CONFIG:Release>:-O3 -DNDEBUG>
  $<$<CONFIG:Profile>:-O3 -ggdb3 -DNDEBUG -fno-inline>
  $<$<CONFIG:Coverage>:-O0 -g --coverage>
  $<IF:$<BOOL:${PORTABLE_BUILD}>,-march=x86-64,-march=native>
)

# No TBB: the code is fully sequential (no std::execution / parallel STL). The earlier
# unconditional TBB link was a vestige of a removed std::execution::par query path (commit
# 7f9b957 made queries fully sequential); linking it only added an unused shared-lib dependency.

add_library(warning_flags INTERFACE)
target_compile_options(warning_flags INTERFACE
  "-Wall"
  "-Wextra"
  "-pedantic"
  $<$<CONFIG:Release>:
   $<$<BOOL:${UNNAMED_WERROR}>:-Werror>>
)

add_library(asan_flags INTERFACE)
target_compile_options(asan_flags INTERFACE
  $<$<CONFIG:Debug>:-fno-omit-frame-pointer -fsanitize=address>
)
target_link_options(asan_flags INTERFACE
  $<$<CONFIG:Debug>:-fno-omit-frame-pointer -fsanitize=address>
)
# Only GCC needs the explicit libasan link; with Clang, -fsanitize=address (above) already pulls in
# compiler-rt's ASan runtime, and additionally linking GCC's libasan triggers an
# "incompatible ASan runtimes" abort. -fsanitize=address handles the runtime for both compilers.
target_link_libraries(asan_flags INTERFACE
  $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:GNU>>:asan>
)

# add_subdirectory(thirdparty)
add_subdirectory(lib)
add_subdirectory(app)
# add_subdirectory(benchmark)
add_subdirectory(thirdparty)

if (WITH_UNION_BENCH)
  add_subdirectory(benchmark/union_bench)
endif()

if (WITH_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
