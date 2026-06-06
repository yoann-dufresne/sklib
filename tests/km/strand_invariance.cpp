// v0.4.2 — strand-invariant per-k-mer canonical framing.
//
// The residual short-minimizer bug: when the minimal central m-mer (the minimizer) repeats inside
// a k-mer, the sliding window picked a non-mirror-symmetric occurrence, so a k-mer stored on one
// strand was framed at a different column when queried on the other -> false negative. The fix
// re-derives each k-mer's frame from its own nucleotides (minimal phi-rank; tie -> most central;
// tie -> smaller full canonical interleaved value), which is strand-invariant and context-free.
//
// These tests (1) pin the interleaved bit layout of the new encoder/decoder against the manipulator
// that actually produces stored records, and (2) assert the discriminating property: every k-mer
// AND its reverse complement is queryable, across small m where the bug was frequent.

#include <gtest/gtest.h>

#include <cstdint>
#include <string>
#include <vector>
#include <random>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/VirtualSkmer.hpp>

namespace {

using kuint = uint64_t;

std::string revcomp(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r)
        switch (c) {
            case 'A': c = 'T'; break;  case 'T': c = 'A'; break;
            case 'C': c = 'G'; break;  case 'G': c = 'C'; break;
        }
    return r;
}

// nucleotides[] in Skmer.hpp is {A,C,T,G} = values {0,1,2,3} via (ascii>>1)&3.
std::vector<kuint> encode(const std::string& s) {
    std::vector<kuint> v;
    for (char c : s) v.push_back(static_cast<kuint>((c >> 1) & 0b11));
    return v;
}

bool kmer_is_queryable(km::sortedlist::SortedVirtualSkmerList<kuint>& list,
                       km::SkmerManipulator<kuint>& manip, std::string kmer) {
    km::SeqSkmerator<kuint> qg{manip, kmer};
    std::vector<km::Skmer<kuint>> query;
    for (const km::Skmer<kuint>& s : qg) query.push_back(s);
    for (const std::vector<uint8_t>& row : list.query_skmer_batch(query))
        for (uint8_t v : row)
            if (v != 0) return true;
    return false;
}

void expect_strand_invariant(uint64_t k, uint64_t m, const std::string& seq) {
    km::SkmerManipulator<kuint> manip{k, m};
    std::string s = seq;
    km::SeqSkmerator<kuint> gen{manip, s};
    std::vector<km::Skmer<kuint>> enumeration;
    for (const km::Skmer<kuint>& sk : gen) enumeration.push_back(sk);

    km::sortedlist::SortedVirtualSkmerList<kuint> list(k, m);
    list.generate_sorted_list_from_enumeration(enumeration);

    for (size_t i = 0; i + k <= seq.size(); i++) {
        const std::string kmer = seq.substr(i, k);
        EXPECT_TRUE(kmer_is_queryable(list, manip, kmer))
            << "FWD k=" << k << " m=" << m << " seq=" << seq << " kmer[" << i << "]=" << kmer;
        EXPECT_TRUE(kmer_is_queryable(list, manip, revcomp(kmer)))
            << "RC  k=" << k << " m=" << m << " seq=" << seq << " kmer[" << i << "]=" << kmer
            << " rc=" << revcomp(kmer);
    }
}

}  // namespace

// (1a) The encoder/decoder bit layout must match the layout the manipulator (add_nucleotide)
// actually produces for stored records. Feed a known sequence through add_nucleotide, then
// decode_to_nucleotides of the forward window must return the last 2k-m nucleotides fed.
TEST(StrandInvariance, DecodeMatchesAddNucleotideLayout) {
    std::mt19937_64 rng(2026);
    const char nt[4] = {'A', 'C', 'T', 'G'};
    for (uint64_t m : {2UL, 3UL, 5UL, 7UL}) {
        for (uint64_t k : {m + 2, 2 * m, 2 * m + 3}) {
            const uint64_t L = 2 * k - m;
            km::SkmerManipulator<kuint> manip{k, m};
            std::uniform_int_distribution<int> d(0, 3);
            for (int rep = 0; rep < 30; rep++) {
                manip.init_skmer();
                std::string seq;
                for (uint64_t i = 0; i < 3 * L; i++) seq += nt[d(rng)];  // overfill the window
                for (char c : seq) manip.add_nucleotide(static_cast<kuint>((c >> 1) & 0b11));

                km::Skmer<kuint> w = manip.m_fwd;  // forward window of the last L nucleotides
                w.m_pref_size = static_cast<uint16_t>(k - m);
                w.m_suff_size = static_cast<uint16_t>(k - m);

                const std::vector<kuint> got = manip.decode_to_nucleotides(w);
                const std::vector<kuint> expected = encode(seq.substr(seq.size() - L, L));
                ASSERT_EQ(got, expected) << "k=" << k << " m=" << m;
            }
        }
    }
}

// (1b) build_skmer_from_nucleotides is the exact inverse of decode_to_nucleotides for every framing
// position j (full k-mer and partial frames).
TEST(StrandInvariance, EncodeDecodeRoundTrip) {
    std::mt19937_64 rng(99);
    for (uint64_t m : {2UL, 3UL, 5UL, 7UL}) {
        for (uint64_t k : {m + 2, 2 * m, 2 * m + 2}) {
            km::SkmerManipulator<kuint> manip{k, m};
            std::uniform_int_distribution<int> d(0, 3);
            for (int rep = 0; rep < 50; rep++) {
                std::vector<kuint> K(k);
                for (auto& x : K) x = static_cast<kuint>(d(rng));
                for (uint64_t j = 0; j <= k - m; j++) {
                    km::Skmer<kuint> p = manip.build_skmer_from_nucleotides(K, 0, k - 1, j);
                    ASSERT_EQ(manip.decode_to_nucleotides(p), K)
                        << "k=" << k << " m=" << m << " j=" << j;
                }
            }
        }
    }
}

// (2) The fix: every k-mer and its reverse complement is queryable, across small m.
TEST(StrandInvariance, RandomSequences) {
    std::mt19937_64 rng(12345);
    const char nt[4] = {'A', 'C', 'T', 'G'};
    for (uint64_t m : {2UL, 3UL, 4UL, 5UL, 6UL, 7UL}) {
        for (uint64_t k : {m + 3, 2 * m, 2 * m + 1}) {
            std::uniform_int_distribution<int> d(0, 3);
            for (int rep = 0; rep < 4; rep++) {
                std::string seq;
                for (uint64_t i = 0; i < 6 * k; i++) seq += nt[d(rng)];
                expect_strand_invariant(k, m, seq);
            }
        }
    }
}

// (2b) Adversarial tandem repeats of NON-palindromic units force the minimizer to repeat within
// k-mers (the exact non-palindrome trigger of the residual bug).
TEST(StrandInvariance, TandemRepeatNonPalindromeMinimizers) {
    auto tile = [](const std::string& unit, size_t len) {
        std::string s;
        while (s.size() < len) s += unit;
        return s.substr(0, len);
    };
    for (uint64_t m : {2UL, 3UL, 4UL, 5UL, 7UL}) {
        for (uint64_t k : {2 * m, 2 * m + 1, 3 * m}) {
            for (const std::string& unit : {std::string("ACG"), std::string("AAC"),
                                            std::string("ATG"), std::string("AAAC"),
                                            std::string("ACGT"), std::string("AACG")}) {
                expect_strand_invariant(k, m, tile(unit, 6 * k));
            }
        }
    }
}
