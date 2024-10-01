set(EXEC_NAME "rsqf")

file(GLOB_RECURSE DBG_SOURCES "*.cpp")
add_executable(${EXEC_NAME} ${DBG_SOURCES})

target_compile_features(${EXEC_NAME} PUBLIC cxx_std_17)
target_link_libraries(${EXEC_NAME} PUBLIC
  $<BUILD_INTERFACE:build_flags>
  $<BUILD_INTERFACE:warning_flags>
  $<BUILD_INTERFACE:asan_flags>
  dbglib
)

if (STATIC_BUILD)
  target_link_libraries(${EXEC_NAME} PUBLIC
    -Wl,--whole-archive
    pthread
    -Wl,--no-whole-archive
  )
else()
  target_link_libraries(${EXEC_NAME} PUBLIC pthread cxxopts)
endif()

if (ipo_supported)
  set_property(TARGET ${EXEC_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

if (STATIC_BUILD AND NOT ${CMAKE_BUILD_TYPE} MATCHES "Debug")
  target_link_options(${EXEC_NAME} PRIVATE -static)
endif()

# add_dependencies(${EXEC_NAME} deps)

install (TARGETS ${EXEC_NAME}
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
         COMPONENT ${DBG_LIB}_Runtime
)