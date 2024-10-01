
set(DBGLIB_TEST_DATA "${CMAKE_SOURCE_DIR}/tests/data")

file(GLOB_RECURSE TEST_SOURCES "*.cpp")
add_executable(dbglib-tests ${TEST_SOURCES})
target_link_libraries(dbglib-tests PRIVATE build_flags dbglib gtest zlib kseqpp pthread)

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
