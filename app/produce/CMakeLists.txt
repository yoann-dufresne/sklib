set(EXEC_NAME "sskm-produce")

add_executable(${EXEC_NAME} produce_skmers.cpp)

target_compile_features(${EXEC_NAME} PUBLIC cxx_std_23)
target_link_libraries(${EXEC_NAME} PUBLIC
  $<BUILD_INTERFACE:build_type_flags>
  $<BUILD_INTERFACE:warning_flags>
  $<BUILD_INTERFACE:asan_flags>
  sklib
)
target_link_libraries(${EXEC_NAME} PUBLIC pthread)

if (ipo_supported)
  # Release only (see app/sortedlists/CMakeLists.txt).
  set_property(TARGET ${EXEC_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()

install (TARGETS ${EXEC_NAME}
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
         COMPONENT ${SKLIB_LIB}_Runtime
)
