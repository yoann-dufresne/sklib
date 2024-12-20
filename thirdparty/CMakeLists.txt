cmake_policy(SET CMP0077 NEW)

set(THIRD_DIR ${CMAKE_SOURCE_DIR}/thirdparty)

set(EXT_INSTALL_DIR ${CMAKE_BINARY_DIR}/ext-install)

include(FetchContent)
set(FETCHCONTENT_BASE_DIR "${CMAKE_BINARY_DIR}/ext-src")

include(GNUInstallDirs)




find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()






if (WITH_TESTS)
  add_subdirectory(gtest)
endif()


# --- Zlib ---


include(FetchContent)

# Download and configure Zlib
FetchContent_Declare(
    zlib
    URL https://github.com/madler/zlib/archive/refs/tags/v1.2.13.tar.gz
)

FetchContent_MakeAvailable(zlib)


include(ExternalProject)

# --- kseq lib ---

add_subdirectory(kseqpp)
add_library(kseqpp INTERFACE)
target_include_directories(kseqpp INTERFACE ${THIRD_DIR}/kseqpp/include ${THIRD_DIR}/kseqpp/include/kseq++)


# --- parse opt ---

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/thirdparty/cxxopts/CMakeLists.txt")
    message(FATAL_ERROR "cxxopts submodule not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

add_library(cxxopts INTERFACE)

ExternalProject_Add(
  CXXOPTS
  PREFIX CXXOPTS
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/thirdparty/cxxopts
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)

target_include_directories(cxxopts INTERFACE ${THIRD_DIR}/cxxopts/include)