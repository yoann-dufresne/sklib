configure_file("${CMAKE_SOURCE_DIR}/lib/config.hpp.in" "${CMAKE_BINARY_DIR}/include/sklib/config.hpp")

set(SKLIB_LIB "sklib")
# if (NOT APPLE)
#   set(CMAKE_INSTALL_RPATH $ORIGIN)
# endif()


file(GLOB_RECURSE SKLIB_SOURCES "src/*.cpp")
add_library(${SKLIB_LIB} ${SKLIB_SOURCES})

set_target_properties(${SKLIB_LIB} PROPERTIES OUTPUT_NAME "sk")
target_compile_definitions(${SKLIB_LIB} PUBLIC SKLIB_PUBLIC)

add_library(${SKLIB_LIB}::${SKLIB_LIB} ALIAS ${SKLIB_LIB})

target_compile_features(${SKLIB_LIB} PUBLIC cxx_std_20)

target_link_libraries(${SKLIB_LIB} PUBLIC
  $<BUILD_INTERFACE:build_type_flags>
  $<BUILD_INTERFACE:warning_flags>
  zlib pthread kseqpp
)

target_include_directories(${SKLIB_LIB} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/lib/include>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set_target_properties(${SKLIB_LIB} PROPERTIES
  SOVERSION ${CMAKE_PROJECT_VERSION_MAJOR}
  VERSION ${CMAKE_PROJECT_VERSION}
)

# LTO (Release only): emit ColinearChaining.cpp as LTO bitcode so colinear_chaining (called per
# column in the phase-2 compaction hotspot) can inline into the apps that link sklib. `ipo_supported`
# is set in the root CMakeLists after probing the toolchain. Release-gated so Debug/Profile stay fast.
if (ipo_supported)
  set_property(TARGET ${SKLIB_LIB} PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()
