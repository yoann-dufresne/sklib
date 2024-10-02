file(GLOB_RECURSE DBGLIB_SOURCES "src/datastruct/*.cpp")
configure_file("${CMAKE_SOURCE_DIR}/lib/config.hpp.in" "${CMAKE_BINARY_DIR}/include/dbglib/config.hpp")
set(DBGLIB_LIB "dbglib")
# if (NOT APPLE)
#   set(CMAKE_INSTALL_RPATH $ORIGIN)
# endif()
add_library(${DBGLIB_LIB} ${DBGLIB_SOURCES})
add_library(${DBGLIB_LIB}::${DBGLIB_LIB} ALIAS ${DBGLIB_LIB})

set_target_properties(${DBGLIB_LIB} PROPERTIES 
  # OUTPUT_NAME "dbg"
  VERSION ${CMAKE_PROJECT_VERSION}
  SOVERSION ${CMAKE_PROJECT_VERSION_MAJOR}
)

# target_compile_features(${DBGLIB_LIB} PUBLIC cxx_std_17)

# target_compile_definitions(${DBGLIB_LIB} PUBLIC DBGLIB_PUBLIC)

target_link_libraries(${DBGLIB_LIB} PRIVATE
  $<BUILD_INTERFACE:build_flags>
  $<BUILD_INTERFACE:warning_flags>
  zlib pthread kseqpp
)

target_include_directories(${DBGLIB_LIB} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/lib/include>
  # $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# add_dependencies(${KMINDEX_LIB} deps)