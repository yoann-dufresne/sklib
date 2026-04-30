#include <iostream>
#include <cstdlib>
#include <gtest/gtest.h>
#include <string>

#include <io/Kmerator.hpp>

using namespace std;


static const std::string data_path {std::getenv("SKLIB_TEST_DATA") == nullptr ? "" : std::getenv("SKLIB_TEST_DATA")};

// TODO: tmp usage of relative path befor solving the getenv issue.
static const std::string fa0 {"../tests/data/fasta0.fa"};


TEST(Kmerator, init)
{
    using kuint = uint64_t;

    km::KmerManipulator<kuint> manip {5};
    km::FileKmerator<kuint> kmerator {fa0, manip};

    uint64_t nb_kmer {0};
    for ([[maybe_unused]]km::Kmer<kuint> kmer : kmerator)
        nb_kmer += 1;
 
    EXPECT_EQ(nb_kmer, 24);
}
