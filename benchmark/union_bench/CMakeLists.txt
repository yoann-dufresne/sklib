# Isolated mono-thread benchmark/verify harness for set_union. Built only when WITH_UNION_BENCH=ON
# (off by default, so normal builds are unaffected). Mirrors app/sortedlists for flags/links.
set(EXEC_NAME "union_bench")

add_executable(${EXEC_NAME} union_bench.cpp)

target_compile_features(${EXEC_NAME} PUBLIC cxx_std_23)
target_link_libraries(${EXEC_NAME} PUBLIC
  $<BUILD_INTERFACE:build_type_flags>
  $<BUILD_INTERFACE:warning_flags>
  sklib
  pthread
)

if (ipo_supported)
  set_property(TARGET ${EXEC_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
