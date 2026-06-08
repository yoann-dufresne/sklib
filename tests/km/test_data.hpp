#ifndef SKLIB_TEST_DATA_HPP
#define SKLIB_TEST_DATA_HPP

#include <string>

// Resolves a fixture file name (e.g. "fasta0.fa") to its full path under tests/data. The directory
// is injected by CMake as SKLIB_TEST_DATA_DIR (see tests/CMakeLists.txt) so the suite finds its
// data regardless of the working directory it is launched from — previously the tests used a
// build-dir-relative "../tests/data/..." path that only resolved when run from inside build/.
#ifndef SKLIB_TEST_DATA_DIR
#define SKLIB_TEST_DATA_DIR "../tests/data"
#endif

inline std::string test_data(const std::string& name) {
    return std::string(SKLIB_TEST_DATA_DIR) + "/" + name;
}

#endif // SKLIB_TEST_DATA_HPP
