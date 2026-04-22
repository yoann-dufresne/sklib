
set(DBGLIB_TEST_DATA "${CMAKE_SOURCE_DIR}/tests/data")

file(GLOB_RECURSE TEST_SOURCES "*.cpp")
# add_library(${DBGLIB_LIB} ${DBGLIB_SOURCES})
add_executable(dbglib-tests ${TEST_SOURCES})

target_link_libraries(dbglib-tests PUBLIC gtest gtest_main zlib kseqpp pthread dbglib)

if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0.0)
    target_link_libraries(dbglib-tests PUBLIC stdc++fs)
endif()


add_test(
    NAME dbglib:all
    COMMAND dbglib-tests
)

# TODO: Why this is not working ?
# set_property(
#   TEST dbglib:all
#   PROPERTY ENVIRONMENT "DBGLIB_TEST_DATA=${CMAKE_SOURCE_DIR}/tests/data"
# )


file(GLOB TEST_SCRIPTS "${CMAKE_SOURCE_DIR}/tests/*.sh")

foreach(script ${TEST_SCRIPTS})
    configure_file(${script} ${CMAKE_BINARY_DIR} COPYONLY)
endforeach()