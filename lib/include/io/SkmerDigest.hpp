#include <cstdint>
#include <cstddef>
#include <io/Skmer.hpp>

#ifndef SKMER_DIGEST_H
#define SKMER_DIGEST_H

namespace km
{

// Order-sensitive FNV-1a rolling fingerprint over a stream of yielded super-k-mers.
//
// It folds ONLY the semantic fields of each record -- the two interleaved pair words and the
// prefix/suffix sizes -- and never the Skmer tail padding (m_pad). The digest is therefore
// independent of the compiler-chosen struct layout and reproducible across builds: two skmer
// streams hash equal iff they are the same sequence of (value, pref, suff) triples in the same
// order. Any divergence (a changed value, a different framing, a reordering, a missing or extra
// record) flips the digest.
//
// This is the bit-exact regression contract used to guard monothread optimizations of the
// producer (SeqSkmerator / SkmerManipulator): the producer's output must stay identical, so a
// single 64-bit digest over the whole stream is enough to certify equivalence before vs after an
// optimization. The isolated `sskm-produce` binary and the `skmerator_digest` test share this one
// definition so their fingerprints are directly comparable and cannot drift apart.
//
// Note: this is a deliberate, brittle-by-design *regression* check; it is the complement of the
// order-/φ-agnostic coverage invariant the other Skmerator tests assert (which survives intended
// algorithm changes). Regenerate the golden values only when the producer's output semantics are
// changed on purpose.

inline constexpr uint64_t SKMER_DIGEST_INIT  = 1469598103934665603ull; // FNV-1a 64-bit offset basis
inline constexpr uint64_t SKMER_DIGEST_PRIME = 1099511628211ull;       // FNV-1a 64-bit prime

inline void skmer_digest_fold_bytes(uint64_t& h, const void* data, std::size_t n)
{
    const unsigned char* b = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < n; ++i)
        h = (h ^ b[i]) * SKMER_DIGEST_PRIME;
}

// Fold one yielded super-k-mer into the running digest `h` (seed it with SKMER_DIGEST_INIT).
template <typename kuint>
inline void skmer_digest_update(uint64_t& h, const km::Skmer<kuint>& sk)
{
    // The two interleaved words hold the whole packed super-k-mer value; fold the raw object
    // bytes (same approach as Skmer::pair_hasher -- well-defined for every kuint backend, _BitInt
    // included). On the x86_64 little-endian targets this library builds for, the byte order is
    // fixed, so the digest is stable across runs and machines.
    const kuint w0 {sk.m_pair.m_value[0]};
    const kuint w1 {sk.m_pair.m_value[1]};
    skmer_digest_fold_bytes(h, &w0, sizeof(kuint));
    skmer_digest_fold_bytes(h, &w1, sizeof(kuint));

    const uint16_t pref {sk.m_pref_size};
    const uint16_t suff {sk.m_suff_size};
    skmer_digest_fold_bytes(h, &pref, sizeof(pref));
    skmer_digest_fold_bytes(h, &suff, sizeof(suff));
}

} // namespace km

#endif // SKMER_DIGEST_H
