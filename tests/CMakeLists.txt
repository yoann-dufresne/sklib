
set(SKLIB_TEST_DATA "${CMAKE_SOURCE_DIR}/tests/data")

file(GLOB_RECURSE TEST_SOURCES "*.cpp")
# add_library(${SKLIB_LIB} ${SKLIB_SOURCES})
add_executable(sklib-tests ${TEST_SOURCES})

target_link_libraries(sklib-tests PUBLIC gtest gtest_main zlib kseqpp pthread sklib)

if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0.0)
    target_link_libraries(sklib-tests PUBLIC stdc++fs)
endif()


add_test(
    NAME sklib:all
    COMMAND sklib-tests
)

# TODO: Why this is not working ?
# set_property(
#   TEST sklib:all
#   PROPERTY ENVIRONMENT "SKLIB_TEST_DATA=${CMAKE_SOURCE_DIR}/tests/data"
# )


file(GLOB TEST_SCRIPTS "${CMAKE_SOURCE_DIR}/tests/*.sh")

foreach(script ${TEST_SCRIPTS})
    configure_file(${script} ${CMAKE_BINARY_DIR} COPYONLY)
endforeach()
