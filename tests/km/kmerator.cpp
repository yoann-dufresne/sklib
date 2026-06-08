#include <iostream>
#include <cstdlib>
#include <gtest/gtest.h>
#include <string>

#include <io/Kmerator.hpp>
#include "test_data.hpp"

using namespace std;


static const std::string fa0 {test_data("fasta0.fa")};


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
