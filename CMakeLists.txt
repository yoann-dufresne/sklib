cmake_minimum_required(VERSION 3.22.0)
project(sklib VERSION 0.5.0 LANGUAGES C CXX)

set(PROJECT_DESCRIPTION "Sorted k-mer library — compact sort and query of k-mers")

# set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

cmake_policy(SET CMP0135 NEW)

option(PORTABLE_BUILD "Generic x86_64 build" OFF)
option(WITH_APP "Compile CLI applications" ON)
option(WITH_TESTS "Compile tests" ON)
option(STATIC_BUILD "Static build" OFF)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE="Release")
endif()

set(CMAKE_CXX_STANDARD 23)

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
target_link_libraries(asan_flags INTERFACE
  $<$<CONFIG:Debug>:asan>
)

# add_subdirectory(thirdparty)
add_subdirectory(lib)
add_subdirectory(app)
# add_subdirectory(benchmark)
add_subdirectory(thirdparty)

if (WITH_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
