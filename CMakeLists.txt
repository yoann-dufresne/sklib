cmake_minimum_required(VERSION 3.22.0)
project(dbglib VERSION 0.0.1 LANGUAGES CXX)

set(PROJECT_DESCRIPTION "dbglib")

# set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

option(DBGLIB_X86_64 "x86-64 build" OFF)
option(DBGLIB_ARM "ARM build" OFF)
option(PORTABLE_BUILD "Generic x86_64 build" OFF) # only works with DBGLIB_X86_64
option(WITH_APP "Compile kmindex cli" ON)
option(WITH_TESTS "Compile tests" ON)
option(STATIC_BUILD "Static build" OFF)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE="Release")
endif()

set(CMAKE_CXX_STANDARD 20)

add_library(build_flags INTERFACE)
target_compile_options(build_flags INTERFACE
  $<$<CONFIG:Debug>:-O -g>
  $<$<CONFIG:Release>:-O3 -DNDEBUG>
  $<$<CONFIG:Profile>:-O3 -ggdb3 -DNDEBUG -fno-inline>
  $<$<CONFIG:Coverage>:-O0 -g --coverage>
  $<$<BOOL:${DBGLIB_X86_64}>: $<IF:$<BOOL:${PORTABLE_BUILD}>,-march=x86-64,-march=native> -mbmi2 -msse4.2> # x86-64-specific options
  $<$<BOOL:${BIOLIB_ARM}>: > # ARM-specific options
)

# target_compile_definitions(FOR_X86_64_PLATFORM)
add_compile_definitions(FOR_ARM_PLATFORM)

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

add_subdirectory(lib)
add_subdirectory(app)
add_subdirectory(thirdparty)

if (WITH_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
