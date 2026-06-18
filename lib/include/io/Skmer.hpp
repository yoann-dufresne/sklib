#include <cstdint>
#include <iostream>
#include <array>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <io/wide_int.hpp>

#ifndef SKMER_H
#define SKMER_H


using namespace std;


namespace km
{


/** This superkmer class represents a full size superkmer (ie 2*k-m nucleotides). It is stored in a interleavec way:
 * It means the the most significant nucleotide in the middle one in the superkmer. For example if you have a superkmer
 * ABCDE, its representation here will be CDBEA. We made this choice as the more a nucleotide is "central" in the
 * superkmer, the most it is shared by the kmers from the superkmers.
 **/
template<typename kuint>
class Skmer
{

public:
    struct pair;

    pair m_pair;
    uint16_t m_pref_size;
    uint16_t m_suff_size;
    // Explicit, zero-initialized tail padding. The serializer dumps records raw, so any byte the
    // compiler leaves as alignment padding leaks uninitialized heap and makes the output bytes (not
    // the k-mer set) differ run-to-run. A single uint32_t only covered the narrow widths: at
    // kuint=__uint128_t the record aligns to 16, the data ends at 40, and the 8 trailing bytes stayed
    // uninitialized (nondeterministic on-disk records, e.g. across construction thread counts). Size
    // the pad from kuint (pair is still incomplete here; sizeof(pair)==2*sizeof(kuint),
    // alignof(pair)==alignof(kuint)) so the named members fill the record to its natural size for
    // EVERY width, with no compiler-inserted tail padding — and sizeof is unchanged (4 bytes at 32/64
    // bit, 12 at 128), so the on-disk layout is identical.
    static constexpr size_t k_data_bytes   = 2 * sizeof(kuint) + 2 * sizeof(uint16_t);
    static constexpr size_t k_record_bytes =
        ((k_data_bytes + sizeof(uint32_t) + alignof(kuint) - 1) / alignof(kuint)) * alignof(kuint);
    std::array<std::uint8_t, k_record_bytes - k_data_bytes> m_pad{};

    Skmer() : m_pair(), m_pref_size(0), m_suff_size(0)
    {}
    Skmer(const pair& value) : m_pair(value), m_pref_size(0), m_suff_size(0)
    {}
    Skmer(const pair& value, uint16_t prefix_size, uint16_t suffix_size) : m_pair(value), m_pref_size(prefix_size), m_suff_size(suffix_size)
    {}
    Skmer(const Skmer& other) : m_pair(other.m_pair), m_pref_size(other.m_pref_size), m_suff_size(other.m_suff_size)
    {}
    Skmer(const Skmer&& other) : m_pair(other.m_pair), m_pref_size(other.m_pref_size), m_suff_size(other.m_suff_size)
    {}

    Skmer<kuint>& operator= (const Skmer<kuint>& other)
    { m_pair = other.m_pair; m_pref_size = other.m_pref_size; m_suff_size = other.m_suff_size; return *this; }

    Skmer<kuint>& operator= (const pair& value)
    { m_pair = value; return *this; }

    Skmer<kuint>& operator= (Skmer<kuint>&& other)
    { m_pair = other.m_pair; m_pref_size = other.m_pref_size; m_suff_size = other.m_suff_size; return *this; }

    bool operator<(const Skmer<kuint>& other) const
    {
        return m_pair < other.m_pair;
    }

    bool operator<=(const Skmer<kuint>& other) const
    {
        return m_pair <= other.m_pair;
    }

    bool operator==(const Skmer<kuint>& other) const
    {
        return (m_pair == other.m_pair && m_pref_size == other.m_pref_size && m_suff_size == other.m_suff_size) ;
    }


    friend std::ostream& operator<<(std::ostream& os, const Skmer<kuint>& p)
    {
        os << p.m_pair << " pref:" << p.m_pref_size << " suff:" << p.m_suff_size;
        return os;
    }

    /** Neasted struct to manage pair of uints
     **/
    struct pair
    {
        // uint 0 is the less significant
        std::array<kuint, 2> m_value;

        pair() : m_value({0, 0})
        {}
        pair(kuint& single) : m_value({single, 0})
        {}
        pair(kuint single) : m_value({single, 0})
        {}
        pair(const kuint& less_significant, const kuint& most_significant)
                                        : m_value({less_significant, most_significant})
        {}
        pair(const pair& other) : m_value({other.m_value[0], other.m_value[1]})
        {}

        pair& operator= (const pair& other)
        {
            m_value[0] = other.m_value[0];
            m_value[1] = other.m_value[1];

            return *this;
        }

        bool operator<(const pair& other) const
        {
            if (m_value[1] == other.m_value[1])
                return m_value[0] < other.m_value[0];
            else
                return m_value[1] < other.m_value[1];
        }

        bool operator==(const pair& other) const
        { return m_value[0] == other.m_value[0] and m_value[1] == other.m_value[1]; }

        bool operator<=(const pair& other) const
        {
            if (m_value[1] == other.m_value[1])
                return m_value[0] <= other.m_value[0];
            else
                return m_value[1] < other.m_value[1];
        }

        // Without this, `a > b` falls back to the implicit `operator uint64_t()`
        // conversion below and compares only the low word (m_value[0]), which made
        // kmer_compare report distinct k-mers as equal once a k-mer spans both
        // words (issue #5: query false positives for 2*(2k-m) > 64).
        bool operator>(const pair& other) const
        {
            if (m_value[1] == other.m_value[1])
                return m_value[0] > other.m_value[0];
            else
                return m_value[1] > other.m_value[1];
        }

        // Three-way lexicographic compare, most-significant word first: -1 / 0 / +1. A single
        // hi-then-lo pass does at most two word compares, where `operator<` followed by `operator>`
        // (as kmer_compare uses) does up to four. Its sign is identical to that ordered pair, so a
        // caller that switches to it stays byte-identical. Used by the wide-store set-op merge inner
        // loop (SetOperations.hpp merge_columns), where each word is a multi-limb _BitInt.
        int compare3(const pair& other) const
        {
            if (m_value[1] != other.m_value[1]) return m_value[1] < other.m_value[1] ? -1 : 1;
            if (m_value[0] != other.m_value[0]) return m_value[0] < other.m_value[0] ? -1 : 1;
            return 0;
        }

        pair operator~ () const
        {
            return pair(~m_value[0], ~m_value[1]);
        }

        pair& operator>>= (const uint64_t shift)
        {
            if (shift == 0) return *this;
            if (shift >= 2 * sizeof(kuint) * 8)
            {
                // The shift clears the whole pair. Handle it explicitly: the branch below would
                // otherwise evaluate `m_value[1] >> (shift - sizeof(kuint)*8)` with a count >= the
                // kuint width, which is undefined behaviour (and GCC and Clang disagree on the
                // result). Saturating to 0 keeps the operator well-defined for any shift.
                m_value[0] = static_cast<kuint>(0);
                m_value[1] = static_cast<kuint>(0);
            }
            else if (shift >= sizeof(kuint) * 8)
            {
                m_value[0] = static_cast<kuint>(m_value[1] >> (shift - sizeof(kuint) * 8));
                m_value[1] = static_cast<kuint>(0);
            }
            else
            {
                // Most significant kuint
                const kuint transfer_mask {static_cast<kuint>((~static_cast<kuint>(0)) >> (sizeof(kuint) * 8 - shift))};
                const kuint transfer_slice {static_cast<kuint>(m_value[1] & transfer_mask)};
                m_value[1] >>= shift;

                // Less significant kuint
                const uint64_t shift_transfered {sizeof(kuint) * 8 - shift};
                m_value[0] >>= shift;
                m_value[0] |= transfer_slice << shift_transfered;
            }

            return *this;
        }

        pair operator>> (const uint64_t shift) const
        {
            pair p {*this};
            p >>= shift;
            return p;
        }

        pair& operator<<= (const uint64_t shift)
        {
            if (shift == 0) return *this;
            if (shift >= 2 * sizeof(kuint) * 8)
            {
                // The shift clears the whole pair. Handle it explicitly: the branch below would
                // otherwise evaluate `m_value[0] << (shift - sizeof(kuint)*8)` with a count >= the
                // kuint width, which is undefined behaviour (and GCC and Clang disagree on the
                // result). Saturating to 0 keeps the operator well-defined for any shift.
                m_value[0] = static_cast<kuint>(0);
                m_value[1] = static_cast<kuint>(0);
            }
            else if (shift >= sizeof(kuint) * 8)
            {
                m_value[1] = static_cast<kuint>(m_value[0] << (shift - sizeof(kuint) * 8));
                m_value[0] = static_cast<kuint>(0);
            }
            else
            {
                // Less significant kuint
                const kuint transfer_mask {static_cast<kuint>((~static_cast<kuint>(0)) << (sizeof(kuint) * 8 - shift))};
                const kuint transfer_slice {static_cast<kuint>(m_value[0] & transfer_mask)};
                m_value[0] <<= shift;

                // Most significant kuint
                const uint64_t shift_transfered {sizeof(kuint) * 8 - shift};
                m_value[1] <<= shift;
                m_value[1] |= transfer_slice >> shift_transfered;
            }
            return *this;
        }

        pair operator<< (const uint64_t shift) const
        {
            pair p {*this};
            p <<= shift;
            return p;
        }

        pair& operator&= (const pair& other)
        {
            m_value[0] &= other.m_value[0];
            m_value[1] &= other.m_value[1];
            return *this;
        }

        pair operator& (const uint64_t value) const
        {
            pair p{*this};
            p.m_value[0] &= value;
            p.m_value[1] = 0;
            return p;
        }

        pair operator& (const pair& other) const
        {
            pair p{*this};
            p.m_value[0] &= other.m_value[0];
            p.m_value[1] &= other.m_value[1];
            return p;
        }

        pair operator| (const pair& other) const
        {
            pair p{*this};
            p.m_value[0] |= other.m_value[0];
            p.m_value[1] |= other.m_value[1];
            return p;
        }

        pair operator^ (const pair& other) const
        {
            pair p{*this};
            p.m_value[0] ^= other.m_value[0];
            p.m_value[1] ^= other.m_value[1];
            return p;
        }

        pair& operator|= (const pair& other)
        {
            m_value[0] |= other.m_value[0];
            m_value[1] |= other.m_value[1];
            return *this;
        }

        pair& operator|= (const kuint& value)
        {
            m_value[0] |= value;
            return *this;
        }

        operator uint64_t()
        {
            return static_cast<uint64_t>(m_value[0]);
        }

        // Full-width low word, WITHOUT the 64-bit narrowing of operator uint64_t(). Use this
        // wherever a value that may span more than 64 bits (e.g. a 2m-bit minimizer with m >= 33)
        // is extracted from the low word: static_cast<kuint>(pair) would otherwise route through
        // operator uint64_t() and silently truncate to 64 bits (the m>=33 single-bucket bug).
        kuint to_kuint() const
        {
            return m_value[0];
        }


        friend std::ostream& operator<<(std::ostream& os, const typename Skmer<kuint>::pair& p)
        {
            // [kuint 1] Prints the bits one at a time from the most significant to the less significant
            for (size_t idx{sizeof(kuint) * 8} ; idx>0 ; idx--)
                os << ((static_cast<uint64_t>(p.m_value[1] >>  (idx - 1))) & 1);

            os << " ";

            // [kuint 2] Prints the bits one at a time from the most significant to the less significant
            for (size_t idx{sizeof(kuint) * 8} ; idx>0 ; idx--)
                os << ((static_cast<uint64_t>(p.m_value[0] >> (idx - 1))) & 1);

            return os;
        }

    };

    struct pair_hasher
    {
        std::size_t operator()(pair const & p) const {
            // Width-agnostic FNV-1a over the raw word bytes. std::hash is not provided for every
            // kuint backend (notably _BitInt), and pair has no padding (m_value is a contiguous
            // kuint[2]), so folding its bytes is well-defined for any width.
            const unsigned char* bytes = reinterpret_cast<const unsigned char*>(p.m_value.data());
            std::size_t h = 1469598103934665603ull; // FNV-1a offset basis
            for (std::size_t i = 0; i < sizeof(p.m_value); ++i)
                h = (h ^ bytes[i]) * 1099511628211ull; // FNV-1a prime
            return h;
        }
    };

};



/** Keep only the low `eff` bits of a `gen`-width interleaved pair, returned as a `store`-width pair.
 * Used both to down-convert stored skmers/queries and to rebuild a quotient manipulator's masks
 * from a full-width one. Requires `store` no wider than `gen` and `eff <= 2*sizeof(store)*8`.
 **/
template<typename gen, typename store>
inline typename Skmer<store>::pair truncate_pair(const typename Skmer<gen>::pair& s, uint64_t eff)
{
    using gpair = typename Skmer<gen>::pair;
    using spair = typename Skmer<store>::pair;
    const uint64_t sw {sizeof(store) * 8};                          // store word width in bits
    const store lo {static_cast<store>(s.m_value[0])};             // bits [0, sw)
    const store hi {static_cast<store>((gpair(s) >> sw).m_value[0])}; // bits [sw, 2*sw)
    spair out {lo, hi};
    spair mask {static_cast<store>(~static_cast<store>(0)), static_cast<store>(~static_cast<store>(0))};
    mask >>= (2 * sw - eff);                                        // low `eff` bits
    out &= mask;
    return out;
}


template<typename kuint>
class SkmerPrettyPrinter
{
public:
    Skmer<kuint> m_skmer;
    uint64_t k;
    uint64_t m;
    uint64_t sk_size;
    uint64_t m_suff_size;
    uint64_t m_pref_size;

    SkmerPrettyPrinter(uint64_t k, uint64_t m): k(k), m(m)
                    , sk_size(2*k-m), m_suff_size(sk_size / 2), m_pref_size((sk_size+1) / 2)
    {};

    SkmerPrettyPrinter& operator<< (Skmer<kuint> const skmer)
    {
        m_skmer = skmer;
        return *this;
    }

};

template<typename kuint>
std::ostream& operator<<(std::ostream& os, const SkmerPrettyPrinter<kuint> pp)
{
    const uint64_t k {pp.k};
    const uint64_t m {pp.m};
    static const char nucleotides[] = {'A', 'C', 'T', 'G'};
    os << "[skmer not interleaved: ";

    // Forward prefix
    for (uint64_t pref_idx{k-m-pp.m_skmer.m_pref_size} ; pref_idx<pp.m_pref_size ; pref_idx++)
    {
        os << nucleotides[((pp.m_skmer).m_pair >> (4 * pref_idx)) & 0b11UL];
    }
    os << " ";
    // Forward suffix
    for (uint64_t suf_idx{pp.m_suff_size} ; suf_idx>(k-m-pp.m_skmer.m_suff_size) ; suf_idx--)
    {
        os << nucleotides[((pp.m_skmer).m_pair >> (4 * suf_idx - 2)) & 0b11UL];
    }

    os << "]";

    return os;
}

using orientation_t = bool;
const bool forward_c {true};
const bool reverse_c {false};

template<typename kuint>
class SkmerManipulator
{
public:
    using kpair = typename km::Skmer<kuint>::pair;

    const uint64_t k;
    const uint64_t m;
    const uint64_t sk_size;

    Skmer<kuint> m_fwd;
    Skmer<kuint> m_rev;

protected:
    // !!! WARNING !!!
    // m_suff_size and m_pref_size are used for internal representation
    // and used for shifting and masking operations. They are not the same
    // as for superkmers: they do not take into account the minimizer size.
    uint64_t m_suff_size;
    uint64_t m_pref_size;

    bool m_current_orientation;

    Skmer<kuint>::pair m_fwd_suffix_buff;
    Skmer<kuint>::pair m_fwd_prefix_buff;
    Skmer<kuint>::pair m_rev_suffix_buff;
    Skmer<kuint>::pair m_rev_prefix_buff;

    // Precomputed (loop-invariant) word index + in-word offset of the two central-nucleotide
    // transfer slots, so add_nucleotide moves the central nucleotide prefix<->suffix with a scalar
    // word op instead of two branchy full-pair shifts per base. A = topmost suffix slot (bit
    // m_suff_size*4-2), B = topmost prefix slot (bit (m_pref_size-1)*4). The moved value is a single
    // 2-bit nucleotide whose bit pair never straddles a word (both slot bits are ≡0/2 mod 4 and the
    // word width is a multiple of 4), so a per-word read/OR reproduces the pair-shift result exactly.
    uint64_t m_cA_word, m_cA_off, m_cB_word, m_cB_off;

    // Number of high φ-minimizer bits not stored (quotiented into the bucket id). 0 for a
    // full-width list; b>0 narrows the meaningful skmer to 2*(2k-m)-b bits, so m_mask (and every
    // mask derived from it via the generators) covers only the retained low bits. See truncate_skmer.
    const uint64_t m_quotient_bits;

    const Skmer<kuint>::pair max_pair_value;
    const kpair m_mask;

    kpair m_minimizer_mask;
    kpair* m_pref_masks;
    kpair* m_suff_masks;
    // Per-(pref/suff size) "absent flank" fill masks: 0b11 over the absent (low) flank slots. Lets
    // mask_absent_nucleotides (a per-yield hot path, also called inside reverse_complement) be a pair
    // of table lookups instead of two per-call loops. Indexed by prefix/suffix size in [0, k-m].
    kpair* m_absent_pref_masks;
    kpair* m_absent_suff_masks;
    // Per-nibble RC patterns (0x3 / 0xC / 0xA in every nibble): isolate each lane's low/high 2-bit
    // half and the per-nucleotide complement, so reverse_complement swaps prefix<->suffix and
    // complements word-parallel instead of lane-by-lane.
    kuint m_rc_lo;
    kuint m_rc_hi;
    kuint m_rc_compl;

    // --- Minimizer-order permutation φ (hash-prospector-class, bijective on 2m bits) ---
    // The minimizer occupies the top 2m bits of m_pair (above bit 4(k-m)). φ relabels
    // that 2m-bit value via a fixed invertible xorshift-multiply mixer so the order on
    // minimizers no longer favors poly-A (raw value 0). Each step is a bijection on
    // Z/2^{2m}: odd-constant multiply (masked), and x ^= x>>s (in range). A final XOR
    // with a nonzero constant breaks the 0-fixpoint the prospector mixers have
    // (lowbias32(0)=0), so φ(0)≠0 — otherwise poly-A would stay minimal and the whole
    // reordering would be a no-op.
    uint64_t m_phi_w;                  // minimizer width in bits = 2m
    kuint m_phi_mask;                  // (1<<2m)-1, or all-ones if 2m == bits(kuint)
    uint64_t m_phi_s1, m_phi_s2, m_phi_s3; // xorshift amounts, scaled to 2m
    kuint m_phi_c1, m_phi_c2;          // odd multiply constants
    kuint m_phi_c1_inv, m_phi_c2_inv;  // their inverses mod 2^bits(kuint)
    kuint m_phi_k;                     // nonzero XOR constant (ensures φ(0)≠0)


    std::vector<kpair > prefix_suffix_mask;
    std::vector<kpair > kmer_masks;
    std::vector<kpair > nucleotide_masks;

    // // The amount of bit shifts needed to reach the 4 most significant bits of a kuint
    // static constexpr uint64_t uints_middle_shift {sizeof(kuint) * 8 - 4};

public:
    // b = quotient bits: the number of high φ-minimizer bits implied by the bucket id and therefore
    // not stored. Default 0 = full-width list (unchanged behaviour). The retained skmer is the low
    // 2*sk_size - b bits, so m_mask is narrowed accordingly and every generator-built mask follows.
    SkmerManipulator(const uint64_t k, const uint64_t m, const uint64_t b = 0)
        : k(k), m(m), sk_size(2*k-m), m_suff_size(sk_size / 2), m_pref_size((sk_size+1) / 2)
        , m_current_orientation(forward_c)
        // Central-nucleotide transfer slots (members m_cA_*/m_cB_*). Computed here in the init list
        // (not the body) because generate_masks_*() below call add_nucleotide, which reads them, while
        // still initializing prefix_suffix_mask/kmer_masks/nucleotide_masks -- so they must be ready
        // before those members. Declared before the mask vectors, so declaration order guarantees it.
        , m_cA_word((m_suff_size * 4 - 2) / (sizeof(kuint) * 8))
        , m_cA_off((m_suff_size * 4 - 2) % (sizeof(kuint) * 8))
        , m_cB_word(((m_pref_size - 1) * 4) / (sizeof(kuint) * 8))
        , m_cB_off(((m_pref_size - 1) * 4) % (sizeof(kuint) * 8))
        , m_quotient_bits(b)
        , max_pair_value(static_cast<kuint>(~static_cast<kuint>(0)), static_cast<kuint>(~static_cast<kuint>(0)))
        , m_mask( max_pair_value >> (2 * sizeof(kuint) * 8 - (2 * sk_size - b)) ), prefix_suffix_mask(generate_masks_sp())
        , kmer_masks(generate_masks_k()), nucleotide_masks(generate_masks_nucleotide())
    {

        // The retained (post-quotient) skmer occupies 2*sk_size - b bits; it must fit the pair.
        assert(2 * sk_size - b <= 2 * sizeof(kuint) * 8);

        // Compute all the possible prefix/suffix masks
        m_pref_masks = new kpair[k-m+1];
        m_suff_masks = new kpair[k-m+1];

        for (uint64_t i{1} ; i<=k-m ; i++)
        {
            // Prefix mask
            static const kpair pref_seed {static_cast<kuint>(0b0011)};
            m_pref_masks[i] = m_pref_masks[i-1] | (pref_seed << (4 * (k - m - i)));

            // Suffix mask
            static const kpair suff_seed {static_cast<kuint>(0b1100)};
            m_suff_masks[i] = m_suff_masks[i-1] | (suff_seed << (4 * (k - m - i)));
        }

        // Precompute the absent-flank fill masks once (same bits mask_absent_nucleotides would set in
        // its loops): for prefix/suffix size p, OR of 0b11 over the low (k-m-p) flank slots.
        m_absent_pref_masks = new kpair[k-m+1];
        m_absent_suff_masks = new kpair[k-m+1];
        for (uint64_t p{0} ; p<=k-m ; p++)
        {
            kpair pm {};
            kpair sm {};
            for (uint64_t i{0} ; i<(k-m-p) ; i++)
            {
                pm |= kpair(static_cast<kuint>(0b11U)) << (4 * i);
                sm |= kpair(static_cast<kuint>(0b11U)) << (4 * i + 2);
            }
            m_absent_pref_masks[p] = pm;
            m_absent_suff_masks[p] = sm;
        }

        // Per-nibble RC patterns: 0x1 in every nibble, scaled to 0x3 / 0xC / 0xA (each fits a nibble,
        // so the multiply never carries across nibbles).
        const kuint ones_nibble {static_cast<kuint>((~static_cast<kuint>(0)) / static_cast<kuint>(0xF))};
        m_rc_lo    = static_cast<kuint>(ones_nibble * static_cast<kuint>(0x3));
        m_rc_hi    = static_cast<kuint>(ones_nibble * static_cast<kuint>(0xC));
        m_rc_compl = static_cast<kuint>(ones_nibble * static_cast<kuint>(0xA));

        // Minimizer mask
        kpair sub_maks = max_pair_value >> (2 * sizeof(kuint) * 8 - 4 * (k - m));
        m_minimizer_mask = m_mask ^ sub_maks;

        // --- φ (minimizer-order permutation) setup ---
        m_phi_w = 2 * m;
        m_phi_mask = (m_phi_w >= sizeof(kuint) * 8)
                        ? static_cast<kuint>(~static_cast<kuint>(0))
                        : static_cast<kuint>((static_cast<kuint>(1) << m_phi_w) - 1);
        // Shifts scaled to the 2m width (a fixed 16 would degenerate to identity for
        // small m). All ≥ 1 so x ^= x>>s stays a bijection.
        m_phi_s1 = std::max<uint64_t>(1, m);
        m_phi_s2 = std::max<uint64_t>(1, m > 1 ? m - 1 : 1);
        m_phi_s3 = std::max<uint64_t>(1, m);
        // splitmix64 finalizer constants (odd). Used full-width; the multiply is masked
        // to 2m bits in phi(), and the inverse mod 2^bits(kuint) masks down to mod 2^2m.
        m_phi_c1 = static_cast<kuint>(0xbf58476d1ce4e5b9ULL);
        m_phi_c2 = static_cast<kuint>(0x94d049bb133111ebULL);
        m_phi_c1_inv = mul_inverse(m_phi_c1);
        m_phi_c2_inv = mul_inverse(m_phi_c2);
        // Golden-ratio odd constant → low bit set → nonzero under any 2m-bit mask.
        m_phi_k = static_cast<kuint>(0x9E3779B97F4A7C15ULL) & m_phi_mask;

        // A quotiented manipulator (b>0) is search-only: the masks generated above used the
        // narrowed m_mask, which corrupts them (add_nucleotide wraps the suffix buffer at the
        // truncated boundary). Rebuild the search masks from a full-width manipulator.
        if (b > 0)
            regenerate_quotient_masks();

        // Skmer and skmer buffers init
        this->init_skmer();
    }

    ~SkmerManipulator()
    {
        delete[] m_pref_masks;
        delete[] m_suff_masks;
        delete[] m_absent_pref_masks;
        delete[] m_absent_suff_masks;
    }

    void init_skmer()
    {
        // Buffers
        m_fwd_suffix_buff = {0, 0};
        m_fwd_prefix_buff = {0, 0};
        m_rev_suffix_buff = {0, 0};
        m_rev_prefix_buff = {0, 0};

        // Skmers
        m_fwd = Skmer<kuint>{};
        m_rev = Skmer<kuint>{};
    }

    bool skmer_equals(const km::Skmer<kuint>& left, const km::Skmer<kuint>& right)
    {
        if ((left.m_pref_size != right.m_pref_size) or (left.m_suff_size != right.m_suff_size))
            return false;

        kpair& pref_mask = m_pref_masks[left.m_pref_size];
        kpair& suff_mask = m_suff_masks[left.m_suff_size];
        kpair pair_mask {pref_mask | m_minimizer_mask | suff_mask};

        kpair left_masked {left.m_pair & pair_mask};
        kpair right_masked {right.m_pair & pair_mask};

        return left_masked == right_masked;
        // return false;
    }

    /** Add a binarized nucleotide (2bits) to the current skmer.
     * @param nucl a 2-bits binarized nucleotide according the current encoding. Warning: The nucleotide is not checked before insertion, so it can destroy the skmer if not correctly formated.
     * @return The current canonical skmer according to the interleaved order.
     **/
    inline Skmer<kuint>& add_nucleotide(const kuint nucl)
    {
        // --- forward prefix ---
        // Shift prefix to the right, then move the central nucleotide (topmost suffix slot A ->
        // topmost prefix slot B). It is a single 2-bit nucleotide that never straddles a word, so a
        // scalar per-word read/OR reproduces the old pair-shift result with no branchy pair shifts.
        // Read the suffix BEFORE it is rewritten below.
        m_fwd_prefix_buff >>= 4;
        const kuint fwd_central_nucl {static_cast<kuint>((m_fwd_suffix_buff.m_value[m_cA_word] >> m_cA_off) & kuint{0b11})};
        m_fwd_prefix_buff.m_value[m_cB_word] |= static_cast<kuint>(fwd_central_nucl << m_cB_off);

        // --- reverse suffix ---
        // Shift the suffix to the right, then move the central nucleotide (topmost prefix slot B ->
        // topmost suffix slot A). Read the prefix BEFORE it is rewritten below.
        m_rev_suffix_buff >>= 4;
        const kuint rev_central_nucl {static_cast<kuint>((m_rev_prefix_buff.m_value[m_cB_word] >> m_cB_off) & kuint{0b11})};
        m_rev_suffix_buff.m_value[m_cA_word] |= static_cast<kuint>(rev_central_nucl << m_cA_off);

        // --- forward suffix ---
        // Shift the suffix
        m_fwd_suffix_buff <<= 4;
        // Add the new nucleotide
        m_fwd_suffix_buff |= nucl << 2;
        // Remove the transfered nucleotide
        m_fwd_suffix_buff &= m_mask;

        // --- reverse prefix ---
        // Shift the prefix
        m_rev_prefix_buff <<= 4;
        // Add the new complement nucleotide
        const kuint compl_nucl{ static_cast<kuint>((nucl + 2) % 4)};

        m_rev_prefix_buff |= compl_nucl;
        // Remove the transfered nucleotide
        m_rev_prefix_buff &= m_mask;

        // --- Merge the interleaved halves ---
        m_fwd = m_fwd_prefix_buff | m_fwd_suffix_buff;
        m_rev = m_rev_prefix_buff | m_rev_suffix_buff;

        if (m_rev < m_fwd)
        {
            m_current_orientation = reverse_c;
            return m_rev;
        }
        else
        {
            m_current_orientation = forward_c;
            return m_fwd;
        }
    }

    inline Skmer<kuint>& add_empty_nucleotide()
    {
        // empty nucleotide
        const kuint nucl = 0b11U;
        // --- forward prefix --- (scalar central-nucleotide transfer A->B; see add_nucleotide)
        m_fwd_prefix_buff >>= 4;
        const kuint fwd_central_nucl {static_cast<kuint>((m_fwd_suffix_buff.m_value[m_cA_word] >> m_cA_off) & kuint{0b11})};
        m_fwd_prefix_buff.m_value[m_cB_word] |= static_cast<kuint>(fwd_central_nucl << m_cB_off);

        // --- reverse suffix --- (scalar central-nucleotide transfer B->A)
        m_rev_suffix_buff >>= 4;
        const kuint rev_central_nucl {static_cast<kuint>((m_rev_prefix_buff.m_value[m_cB_word] >> m_cB_off) & kuint{0b11})};
        m_rev_suffix_buff.m_value[m_cA_word] |= static_cast<kuint>(rev_central_nucl << m_cA_off);

        // --- forward suffix ---
        // Shift the suffix
        m_fwd_suffix_buff <<= 4;
        // Add the new nucleotide
        m_fwd_suffix_buff |= nucl << 2;
        // Remove the transfered nucleotide
        m_fwd_suffix_buff &= m_mask;

        // --- reverse prefix ---
        // Shift the prefix
        m_rev_prefix_buff <<= 4;
        // Add the new empty nucleotide
        m_rev_prefix_buff |= nucl;
        // Remove the transfered nucleotide
        m_rev_prefix_buff &= m_mask;

        // --- Merge the interleaved halves ---
        m_fwd = m_fwd_prefix_buff | m_fwd_suffix_buff;
        m_rev = m_rev_prefix_buff | m_rev_suffix_buff;

        if (m_rev < m_fwd)
        {
            m_current_orientation = reverse_c;
            return m_rev;
        }
        else
        {
            m_current_orientation = forward_c;
            return m_fwd;
        }
    }

    bool is_forward() const
    {
        return m_current_orientation;
    }

    kuint minimizer() const
    {
        return std::min(
            (m_fwd.m_pair >> (4*(k-m))).to_kuint(),
            (m_rev.m_pair >> (4*(k-m))).to_kuint()
        );
    }

    kuint minimizer(const Skmer<kuint>& skmer) const
    {
        return (skmer.m_pair >> (4*(k-m))).to_kuint();
    }

    // Order key of a skmer's (raw) minimizer: φ of the raw-canonical m-mer. The window
    // minimizer is selected by argmin of this rank instead of the raw value, so poly-A
    // (raw 0 → φ(0)=K) is no longer systematically chosen. φ is bijective, so equal
    // ranks ⇔ equal minimizers (the selection's "same minimizer" test stays correct).
    kuint minimizer_rank(const Skmer<kuint>& skmer) const
    {
        return phi(minimizer(skmer));
    }

    // Multiplicative inverse of an odd constant mod 2^bits(kuint) (Newton iteration:
    // each step doubles the number of correct low bits; 6 steps reach ≥ 64 bits).
    static kuint mul_inverse(kuint a)
    {
        kuint x = a; // correct to 3 low bits for odd a
        for (int i = 0; i < 6; i++)
            x *= static_cast<kuint>(2) - a * x;
        return x;
    }

    // φ : fixed invertible permutation of the 2m-bit minimizer value (hash-prospector
    // style xorshift-multiply, masked to 2m bits, + XOR constant so φ(0)≠0).
    kuint phi(kuint x) const
    {
        const kuint mask = m_phi_mask;
        // The mixer only needs the low 2m bits. When 2m <= 64 (every k up to 64, e.g. k=63 -> 2m=62)
        // do it in uint64_t even if kuint is wider (__uint128/_BitInt): each step's surviving bits are
        // the low 2m, and a 64-bit multiply yields exactly the low 64 product bits -- identical after
        // the mask -- while replacing one full 128/256-bit multiply per call (phi runs once per base).
        // Result (and the digest) is unchanged. Wider 2m keeps full-width kuint arithmetic.
        if (sizeof(kuint) > 8 && m_phi_w <= 64)
        {
            const uint64_t m64 {static_cast<uint64_t>(mask)};
            uint64_t u {static_cast<uint64_t>(x) & m64};
            u ^= u >> m_phi_s1;
            u = (u * static_cast<uint64_t>(m_phi_c1)) & m64;
            u ^= u >> m_phi_s2;
            u = (u * static_cast<uint64_t>(m_phi_c2)) & m64;
            u ^= u >> m_phi_s3;
            u ^= static_cast<uint64_t>(m_phi_k);
            return static_cast<kuint>(u & m64);
        }
        x &= mask;
        x ^= x >> m_phi_s1;
        x = (x * m_phi_c1) & mask;
        x ^= x >> m_phi_s2;
        x = (x * m_phi_c2) & mask;
        x ^= x >> m_phi_s3;
        x ^= m_phi_k;
        return x & mask;
    }

    // Inverse of `y = x ^ (x >> s)` on the 2m-bit space, by iterated doubling of s.
    kuint inv_xorshift_right(kuint x, uint64_t s) const
    {
        for (uint64_t sh = s; sh < m_phi_w; sh <<= 1)
            x ^= x >> sh;
        return x;
    }

    // φ⁻¹ : steps of φ in reverse (cold path: ASCII / decoding only).
    kuint phi_inv(kuint x) const
    {
        const kuint mask = m_phi_mask;
        x &= mask;
        x ^= m_phi_k;
        x = inv_xorshift_right(x, m_phi_s3);
        x = (x * m_phi_c2_inv) & mask;
        x = inv_xorshift_right(x, m_phi_s2);
        x = (x * m_phi_c1_inv) & mask;
        x = inv_xorshift_right(x, m_phi_s1);
        return x & mask;
    }

    // Branch-free reversal of all 64 bits of a uint64_t (classic parallel swap network, ~12 ops).
    static uint64_t reverse_bits64(uint64_t x)
    {
        x = ((x >> 1)  & 0x5555555555555555ULL) | ((x & 0x5555555555555555ULL) << 1);
        x = ((x >> 2)  & 0x3333333333333333ULL) | ((x & 0x3333333333333333ULL) << 2);
        x = ((x >> 4)  & 0x0F0F0F0F0F0F0F0FULL) | ((x & 0x0F0F0F0F0F0F0F0FULL) << 4);
        x = ((x >> 8)  & 0x00FF00FF00FF00FFULL) | ((x & 0x00FF00FF00FF00FFULL) << 8);
        x = ((x >> 16) & 0x0000FFFF0000FFFFULL) | ((x & 0x0000FFFF0000FFFFULL) << 16);
        x = (x >> 32) | (x << 32);
        return x;
    }

    // Bit-reversal of the low 2m bits (bit i ↔ bit 2m-1-i), right-aligned in the 2m window.
    // Its own inverse on the 2m-bit space. Composed with φ to spread the *low* (uniform) bits of
    // φ(min) into the high-order positions the bucketing reads: φ(min) is a window-minimum, so its
    // high bits are biased to 0 and bucketing on them under-fills the buckets; reversing puts the
    // uniform low bits on top, restoring a balanced AND still sort-prefix-contiguous bucketing.
    // Runs once per finalized skmer (permute) / per query (routing), never per base.
    // 2m<=64 (every backend up to k~=64, i.e. the whole common range): one 64-bit swap-network
    // reversal then a right-shift to align the window — O(1), replacing the old O(2m) bit loop that
    // dominated phase-1 construction. 2m>64 (k>~64, kuint256) keeps the simple O(2m)<=126 loop; that
    // path is rare and its construction is already a couple of seconds. Output is bit-identical to
    // the loop either way, so stored records (and on-disk indexes) are unchanged.
    kuint reverse_2m(kuint x) const
    {
        const uint64_t w {m_phi_w};
        if (w <= 64)
        {
            // x is φ(min), masked to the low w<=64 bits, so its value lives entirely in the low word.
            const uint64_t u {reverse_bits64(static_cast<uint64_t>(x)) >> (64 - w)};
            return static_cast<kuint>(u);
        }
        kuint r {static_cast<kuint>(0)};
        for (uint64_t i {0}; i < w; ++i)
            r |= static_cast<kuint>((x >> i) & static_cast<kuint>(1)) << (w - 1 - i);
        return r;
    }

    // Replace the raw minimizer in the top 2m bits of m_pair by ψ(minimizer) = reverse_2m(φ(min))
    // (flanks untouched), so a plain m_pair compare orders the list by ψ-minimizer and the bucket
    // prefix (top b bits) carries φ(min)'s uniform low bits. The minimizer SELECTION stays φ-ranked
    // (minimizer_rank, unchanged); only the STORED representation is ψ. No site compares the stored
    // slot against a bare φ value (sort = slot vs slot; bucket = high bits of slot; the only φ_inv
    // consumer is unpermute below), so the φ→ψ change is self-contained.
    void permute_minimizer_slot(Skmer<kuint>& skmer) const
    {
        const kuint raw {(skmer.m_pair >> (4 * (k - m))).to_kuint()};
        const kpair pm {reverse_2m(phi(raw))};
        skmer.m_pair = (skmer.m_pair & (~m_minimizer_mask)) | (pm << (4 * (k - m)));
    }

    // Inverse of permute_minimizer_slot: recover the raw nucleotide minimizer (decode).
    // ψ⁻¹ = φ⁻¹ ∘ reverse_2m (reverse_2m is self-inverse).
    void unpermute_minimizer_slot(Skmer<kuint>& skmer) const
    {
        const kuint pm {(skmer.m_pair >> (4 * (k - m))).to_kuint()};
        const kpair raw {phi_inv(reverse_2m(pm))};
        skmer.m_pair = (skmer.m_pair & (~m_minimizer_mask)) | (raw << (4 * (k - m)));
    }

    /** Replace the bits of absent nucleotides (outside of registered prefix/suffix) by 0b11.
     * @param skmer Super kmer to modify
     **/
    void mask_absent_nucleotides(Skmer<kuint>& skmer) const
    {
        // O(1) precomputed equivalent of the old per-call prefix/suffix fill loops (set 0b11 over the
        // absent flank slots). Identical bits set; pref/suff sizes are always in [0, k-m].
        skmer.m_pair |= m_absent_pref_masks[skmer.m_pref_size];
        skmer.m_pair |= m_absent_suff_masks[skmer.m_suff_size];
    }

    /** Reverse-complement of a super-k-mer in the minimizer-centered interleaved
     * encoding. Position p (from the start) maps to position (2k-m-1-p), so prefix
     * slot p and suffix slot p swap, and every nucleotide is complemented (XOR 0b10).
     * For odd 2k-m the central prefix slot maps to itself (complemented). pref/suff
     * sizes swap. Absent slots are re-masked to 0b11 afterwards.
     * @param skmer Super-k-mer to reverse-complement
     * @return the reverse-complemented super-k-mer
     **/
    Skmer<kuint> reverse_complement(const Skmer<kuint>& skmer) const
    {
        using kpair = typename Skmer<kuint>::pair;
        constexpr uint64_t W {sizeof(kuint) * 8};

        // RC is a within-lane prefix<->suffix swap plus a per-nucleotide complement: the interleaving
        // is built so linear positions p and 2k-m-1-p share a lane (prefix half <-> suffix half), so
        // reverse-complementing is "swap the two 2-bit halves of every nibble, XOR each with 0b10".
        // Do it word-parallel (a few ops/word) instead of looping per lane -- this RC runs once per
        // yield via canonicalize and was the top producer hot spot.
        auto rc_word = [&](kuint w) -> kuint {
            const kuint swapped {static_cast<kuint>(((w & m_rc_lo) << 2) | ((w & m_rc_hi) >> 2))};
            return static_cast<kuint>(swapped ^ m_rc_compl);
        };
        kuint s0 {rc_word(skmer.m_pair.m_value[0])};
        kuint s1 {rc_word(skmer.m_pair.m_value[1])};

        // Odd 2k-m: the central prefix lane maps to itself (complemented) and has no suffix partner.
        // The generic swap pairs it with its absent suffix half, so rewrite that one nibble: keep only
        // compl(central prefix) in the prefix half, drop the suffix half (re-masked away below).
        if (((sk_size + 1) / 2) > (sk_size / 2))
        {
            const uint64_t cbit {4 * (sk_size / 2)};                          // central prefix lane bit
            const kuint in_c {(cbit < W)
                ? static_cast<kuint>((skmer.m_pair.m_value[0] >> cbit)       & kuint{0b11})
                : static_cast<kuint>((skmer.m_pair.m_value[1] >> (cbit - W)) & kuint{0b11})};
            const kuint out_c {static_cast<kuint>(in_c ^ kuint{0b10})};
            if (cbit < W) s0 = static_cast<kuint>((s0 & ~(static_cast<kuint>(0xF) << cbit))       | (out_c << cbit));
            else          s1 = static_cast<kuint>((s1 & ~(static_cast<kuint>(0xF) << (cbit - W))) | (out_c << (cbit - W)));
        }

        kpair out {s0, s1};
        out &= m_mask; // clear lanes beyond the skmer (the per-lane loop left this tail zero)
        Skmer<kuint> rc {out, skmer.m_suff_size, skmer.m_pref_size};
        mask_absent_nucleotides(rc);
        return rc;
    }

    /** Rewrite a (masked) super-k-mer in its canonical orientation: the one whose
     * minimizer-centered interleaved encoding is the smaller. Idempotent, and
     * independent of the surrounding context, so the same k-mer is stored and
     * queried identically (issue #7).
     * @param skmer Super-k-mer to canonicalize in place
     **/
    void canonicalize(Skmer<kuint>& skmer) const
    {
        Skmer<kuint> rc {reverse_complement(skmer)};
        if (rc.m_pair < skmer.m_pair)
            skmer = rc;
    }

    /** Finalize a yielded super-k-mer for the φ-ordered sorted list: pick the raw
     * canonical orientation (reverse_complement works on the raw nucleotide layout),
     * then permute the minimizer slot in place so the stored value is [φ(min) | flanks].
     * From here on the skmer is in "permuted" space and a plain m_pair compare orders by
     * φ(minimizer). Construction and query both reach this via the Skmerator yield, so
     * they stay consistent. Decode with unpermute_minimizer_slot (φ⁻¹).
     * @param skmer Super-k-mer to canonicalize then permute in place
     **/
    void canonicalize_for_sort(Skmer<kuint>& skmer) const
    {
        canonicalize(skmer);
        permute_minimizer_slot(skmer);
    }

    /** Reverse-complement of a minimizer m-mer (m nucleotides, 2 bits each, raw value).
     * Reverses the nucleotide order and complements each (A<->T, C<->G via XOR 0b10).
     **/
    kuint rc_mmer(kuint mmer) const
    {
        kuint r {0};
        for (uint64_t i{0} ; i<m ; i++)
        {
            r = (r << 2) | ((mmer & kuint{0b11}) ^ kuint{0b10});
            mmer >>= 2;
        }
        return r;
    }

    /** A minimizer that equals its own reverse-complement (only possible for even m,
     * e.g. GC, AT, GCGC...). Such a minimizer does NOT fix the super-k-mer's orientation,
     * so the k-mers sharing it need not all be canonical in the whole-super-k-mer
     * orientation — the Skmerator must then split the super-k-mer per k-mer (issue #7
     * completion). For non-palindrome minimizers the minimizer's strand pins a single
     * consistent orientation for every contained k-mer, so no split is needed.
     **/
    bool minimizer_is_rc_palindrome(kuint mmer) const
    {
        return rc_mmer(mmer) == mmer;
    }


    // ===== Strand-invariant per-k-mer canonical framing (v0.4.2) =====
    // Fixes the residual short-minimizer false negatives: when the minimal central m-mer (the
    // minimizer) repeats inside a k-mer, the sliding window's leftmost-style occurrence choice is
    // not mirror-symmetric (the forward and reverse-complement scans pick mirror occurrences ->
    // different storage column -> false negative). The routines below re-derive each k-mer's frame
    // from its own nucleotides using a symmetric rule (minimal phi-rank; tie -> most central; tie
    // -> smaller full canonical interleaved value), so construction-in-context and isolated-query
    // produce identical, strand-invariant records.

    /** Bit shift of linear super-k-mer position `p` (0 = sequence start) inside the interleaved
     * m_pair. The encoding interleaves positions outward from the centre: position `p` sits in
     * prefix slot `p` (bits [4p:4p+1]) for the first ceil((2k-m)/2) positions, otherwise in suffix
     * slot (2k-m-1-p) (bits [4*(2k-m-1-p)+2 : +3]). Decode and encode share this map, so a
     * decode -> encode round-trip is correct by construction (see SkmerManipulator tests). **/
    uint64_t slot_shift(uint64_t p) const
    {
        const uint64_t pref_slots {(sk_size + 1) / 2};                 // ceil((2k-m)/2)
        return (p < pref_slots) ? (4 * p) : (4 * (sk_size - 1 - p) + 2);
    }

    /** Decode a super-k-mer to its present nucleotides in linear 5'->3' order (this orientation):
     * prefix, then minimizer, then suffix. Length L = pref+m+suff; the minimizer is S[pref..pref+m-1]. **/
    std::vector<kuint> decode_to_nucleotides(const Skmer<kuint>& sk) const
    {
        const uint64_t pref {sk.m_pref_size}, suff {sk.m_suff_size};
        const uint64_t L {pref + m + suff};
        const uint64_t offset {k - m - pref};                          // absolute linear start
        std::vector<kuint> S(L);
        for (uint64_t i {0}; i < L; i++)
            S[i] = static_cast<kuint>(sk.m_pair >> slot_shift(offset + i)) & kuint{0b11};
        return S;
    }

    /** Build the forward (this-orientation) super-k-mer covering S[start..end] with the minimizer
     * m-mer at absolute index `mini_abs` (prefix = S[start..mini_abs-1], suffix = S[mini_abs+m..end]).
     * Absent slots are set to 0b11. Result is not yet canonical / phi-permuted. **/
    Skmer<kuint> build_skmer_from_nucleotides(const std::vector<kuint>& S,
                                              uint64_t start, uint64_t end, uint64_t mini_abs) const
    {
        const uint16_t pref {static_cast<uint16_t>(mini_abs - start)};
        const uint16_t suff {static_cast<uint16_t>(end - mini_abs - m + 1)};
        const uint64_t base {k - m - pref};                            // piece's absolute linear start
        const uint64_t L {static_cast<uint64_t>(pref) + m + suff};
        // Pack the interleaved pair word-by-word in two branch-free regimes (prefix half: bit 4p;
        // suffix half: bit 4*(2k-m-1-p)+2), each a single 2-bit nucleotide that never straddles a word
        // -- a scalar per-word OR instead of slot_shift + a full-pair shift per nucleotide. Identical
        // bits; this is the inner loop of the ambiguous-path k-mer reframing (hot at k=63). See decode.
        constexpr uint64_t W {sizeof(kuint) * 8};
        const uint64_t pref_slots {(sk_size + 1) / 2};
        uint64_t split {(pref_slots > base) ? (pref_slots - base) : uint64_t{0}};
        if (split > L) split = L;
        kuint w0 {0};
        kuint w1 {0};
        for (uint64_t i {0}; i < split; i++)
        {
            const uint64_t b {4 * (base + i)};
            const kuint v {static_cast<kuint>(S[start + i])};
            if (b < W) w0 |= static_cast<kuint>(v << b); else w1 |= static_cast<kuint>(v << (b - W));
        }
        for (uint64_t i {split}; i < L; i++)
        {
            const uint64_t b {4 * (sk_size - 1 - (base + i)) + 2};
            const kuint v {static_cast<kuint>(S[start + i])};
            if (b < W) w0 |= static_cast<kuint>(v << b); else w1 |= static_cast<kuint>(v << (b - W));
        }
        Skmer<kuint> out {kpair{w0, w1}, pref, suff};
        mask_absent_nucleotides(out);
        return out;
    }

    /** Minimizer repeat / palindrome test over the decoded nucleotides S[0..L-1] (minimizer at index
     * `pref`): true if the central m-mer equals its own reverse complement, or it occurs (in either
     * orientation) at >=2 positions. Templated on the roll word type so the per-yield roll runs in
     * uint64_t when the m-mer fits (2m<=64, e.g. k=63 -> 2m=62) even if kuint is wider -- then center /
     * rc_center / fwd and every roll compare+shift are 64-bit, not __uint128, for an identical result
     * (all values are <=2m bits). The set {canon, rc_canon} the old code tested equals {center,
     * rc_center} since rc is an involution, so we compare the rolled m-mer against those two directly. **/
    template<typename mword>
    bool mmer_repeats(const uint8_t* S, uint64_t L, uint64_t pref) const
    {
        mword center {0};
        for (uint64_t t {0}; t < m; t++) center |= static_cast<mword>(S[pref + t]) << (2 * t);
        mword rc_center {0};
        {
            mword c {center};
            for (uint64_t t {0}; t < m; t++) { rc_center = (rc_center << 2) | ((c & mword{0b11}) ^ mword{0b10}); c >>= 2; }
        }
        if (center == rc_center) return true;                          // palindrome -> orientation ambiguous
        const uint64_t hi {2 * (m - 1)};
        mword fwd {0};
        for (uint64_t t {0}; t < m; t++) fwd |= static_cast<mword>(S[t]) << (2 * t);
        uint64_t count {0};
        for (uint64_t pos {0}; ; pos++)
        {
            if ((fwd == center || fwd == rc_center) && ++count >= 2) return true;
            if (pos + m >= L) break;
            fwd = (fwd >> 2) | (static_cast<mword>(S[pos + m]) << hi);
        }
        return false;
    }

    /** Cheap, allocation-free ambiguity test on a (raw-canonical) super-k-mer: true if the minimizer
     * is an RC-palindrome (orientation not pinned), or its canonical m-mer value occurs at >=2 of the
     * super-k-mer's m-mer positions (occurrence not pinned). The unique non-palindrome case (false)
     * is the fast path: the whole super-k-mer is already each k-mer's canonical frame. **/
    bool minimizer_is_ambiguous(const Skmer<kuint>& sk) const
    {
        const uint64_t pref {sk.m_pref_size};
        const uint64_t L {pref + m + sk.m_suff_size};
        if (L < m) return false;
        const uint64_t offset {k - m - pref};
        // Decode the present nucleotides once into a stack buffer, reading each from the two pair
        // words with a plain shift (the pair `>>` operator is branch-heavy and was the per-yield
        // hot spot). 2k-m <= pair-bit-width/2, so the buffer and all bit offsets fit.
        uint8_t S[2 * sizeof(kuint) * 8];
        const kuint w_lo {sk.m_pair.m_value[0]};
        const kuint w_hi {sk.m_pair.m_value[1]};
        constexpr uint64_t W {sizeof(kuint) * 8};
        // Decode in two branch-free regimes instead of a per-nucleotide slot_shift: the leading
        // positions live in the prefix half (bit 4p, increasing), the rest in the suffix half
        // (bit 4*(2k-m-1-p)+2). boundary = first linear index whose position reaches pref_slots.
        const uint64_t pref_slots {(sk_size + 1) / 2};
        uint64_t boundary {(pref_slots > offset) ? (pref_slots - offset) : uint64_t{0}};
        if (boundary > L) boundary = L;
        for (uint64_t i {0}; i < boundary; i++)
        {
            const uint64_t b {4 * (offset + i)};
            S[i] = static_cast<uint8_t>((b < W ? (w_lo >> b) : (w_hi >> (b - W))) & kuint{0b11});
        }
        for (uint64_t i {boundary}; i < L; i++)
        {
            const uint64_t b {4 * (sk_size - 1 - (offset + i)) + 2};
            S[i] = static_cast<uint8_t>((b < W ? (w_lo >> b) : (w_hi >> (b - W))) & kuint{0b11});
        }
        // Roll the central m-mer over the decoded buffer to detect a palindrome or a repeat. Run it in
        // uint64_t when the m-mer fits (2m<=64) even if kuint is wider (k=63): every roll compare+shift
        // is then 64-bit instead of __uint128, for an identical boolean. k<=31 keeps kuint (uint64 is
        // already its width, so unchanged). See mmer_repeats.
        if constexpr (sizeof(kuint) > 8)
        {
            if (2 * m <= 64) return mmer_repeats<uint64_t>(S, L, pref);
        }
        return mmer_repeats<kuint>(S, L, pref);
    }

    /** Canonical-minimizer rank of the m-mer at absolute position `p` in S: phi of the smaller of the
     * m-mer's two interleaved orientations. The minimizer occupies the top interleaved bits, so the
     * canonical orientation's minimizer (and thus this rank) is decided by the m-mer alone, independent
     * of any surrounding k-mer window's flanks -- identical to the rank the full-framing path computes
     * (phi(minimizer(canonicalize(framing)))), but a function of the m-mer only, so it memoizes per
     * position. Uses a flank-less framing (prefix=suffix=0) so only the central m nucleotides matter. **/
    kuint mmer_canonical_rank(const std::vector<kuint>& S, uint64_t p) const
    {
        Skmer<kuint> f {build_skmer_from_nucleotides(S, p, p + m - 1, p)};
        canonicalize(f);
        return phi(minimizer(f));
    }

    /** Strand-invariant minimizer occurrence for the k-mer S[start..start+k-1]: (1) minimal
     * phi(canonical m-mer) rank; (2) tie -> most central position (|2j-(k-m)| minimal, mirror-
     * symmetric); (3) tie (mirror pair) -> smaller full canonical interleaved value. Writes the
     * absolute occurrence index and the canonical orientation of that frame. `rank_table[p]` is the
     * memoized per-position rank (mmer_canonical_rank); only the minimal-rank candidates need a full
     * build+canonicalize (a non-minimal position can never win the selection), so the per-(k-mer,
     * candidate) full builds drop from k-m+1 to the few repeated-minimizer occurrences. **/
    void choose_kmer_minimizer(const std::vector<kuint>& S, uint64_t start,
                               const std::vector<kuint>& rank_table,
                               uint64_t& occ_abs, orientation_t& orient) const
    {
        const uint64_t span {k - m};                                   // last minimizer position in the k-mer
        kuint best_rank {rank_table[start]};
        for (uint64_t j {1}; j <= span; j++) best_rank = std::min(best_rank, rank_table[start + j]);
        // (2) most central among minimal-rank positions, then (3) smaller full canonical value. Build
        // the full framing only for minimal-rank positions; the canonical value (with flanks) and the
        // orientation are window-dependent, so those still need the build -- but only here, rarely.
        uint64_t best_j {span + 1};
        uint64_t best_central {~uint64_t{0}};
        Skmer<kuint> best_canon {};
        Skmer<kuint> best_fwd {};
        for (uint64_t j {0}; j <= span; j++)
        {
            if (rank_table[start + j] != best_rank) continue;
            const int64_t sd {static_cast<int64_t>(2 * j) - static_cast<int64_t>(span)};
            const uint64_t d {static_cast<uint64_t>(sd < 0 ? -sd : sd)};
            Skmer<kuint> f {build_skmer_from_nucleotides(S, start, start + k - 1, start + j)};
            Skmer<kuint> c {f};
            canonicalize(c);
            if (best_j > span || d < best_central
                || (d == best_central && c.m_pair < best_canon.m_pair))
            {
                best_central = d;
                best_j = j;
                best_canon = c;
                best_fwd = f;
            }
        }
        occ_abs = start + best_j;
        orient = (reverse_complement(best_fwd).m_pair < best_fwd.m_pair) ? reverse_c : forward_c;
    }

    /** Re-derive the strand-invariant canonical frame of every k-mer in `sk`, then regroup
     * consecutive k-mers that share the same minimizer occurrence and orientation into compact
     * super-k-mers. Each emitted piece is canonical + phi-permuted (ready for the sorted list).
     * Used only on the rare ambiguous super-k-mer; the common case takes the fast path. **/
    void canonical_pieces(const Skmer<kuint>& sk, std::vector<Skmer<kuint>>& out) const
    {
        const std::vector<kuint> S {decode_to_nucleotides(sk)};
        const uint64_t L {S.size()};
        if (L < k) return;
        const uint64_t nk {L - k + 1};
        // Memoize the per-position canonical-minimizer rank once (it depends only on the m-mer at each
        // position, not on the k-mer window), so choose_kmer_minimizer reads it instead of rebuilding +
        // canonicalizing a full framing for every (k-mer, candidate) pair.
        const uint64_t npos {L - m + 1};
        std::vector<kuint> rank_table(npos);
        for (uint64_t p {0}; p < npos; p++) rank_table[p] = mmer_canonical_rank(S, p);
        std::vector<uint64_t> occ(nk);
        std::vector<uint8_t> orient(nk);                               // not vector<bool> (proxy refs)
        for (uint64_t i {0}; i < nk; i++)
        {
            orientation_t o;
            choose_kmer_minimizer(S, i, rank_table, occ[i], o);
            orient[i] = static_cast<uint8_t>(o);
        }
        uint64_t g0 {0};
        for (uint64_t i {1}; i <= nk; i++)
        {
            const bool same {i < nk && occ[i] == occ[g0] && orient[i] == orient[g0]};
            if (!same)
            {
                // group [g0, i-1]: nucleotides S[g0 .. (i-1)+k-1], shared minimizer at occ[g0].
                // Put the piece in the orientation all its k-mers agreed on (orient[g0]) rather than
                // re-canonicalizing the whole piece: with a palindrome (flank-decided) minimizer the
                // whole-piece orientation can differ from the per-k-mer ones, which would store a
                // k-mer in an orientation its isolated query never reproduces.
                Skmer<kuint> piece {build_skmer_from_nucleotides(S, g0, (i - 1) + k - 1, occ[g0])};
                if (orient[g0] == static_cast<uint8_t>(reverse_c))
                    piece = reverse_complement(piece);
                permute_minimizer_slot(piece);
                out.push_back(piece);
                g0 = i;
            }
        }
    }


    /** Compare 2 kmers included in 2 skmers.
     * @param first_skmer First kmer is included in this skmer
     * @param first_kmer_pos Position of the fist kmer in the first skmer. 0 is the kmer that contains the whole prefix.
     * @param second_skmer Second kmer is included in this skmer
     * @param second_kmer_pos Position of the second kmer in the second skmer. 0 is the kmer that contains the whole prefix.
     * @return true if the first kmer is less than the second one
     **/
    bool inline kmer_lt_kmer(const Skmer<kuint>& first_skmer, const uint64_t first_kmer_pos, const Skmer<kuint>& second_skmer, const uint64_t second_kmer_pos) const
    {
        // 1 - Compute the shift needed to only keep informative nucleotides
        uint64_t const first_missing_nucl {std::max(2*first_kmer_pos-1, 2*(k-m-first_kmer_pos))};
        uint64_t const second_missing_nucl {std::max(2*second_kmer_pos-1, 2*(k-m-second_kmer_pos))};
        uint64_t const mask_size {std::max(first_missing_nucl, second_missing_nucl)};

        // 2 - Compare masked kmers
        // check if first skmer shifted right of mask is < second skmer shifted right of mask
        const auto first_kmer {first_skmer.m_pair >> (2 * mask_size)};
        const auto second_kmer {second_skmer.m_pair >> (2 * mask_size)};

        if (first_kmer != second_kmer){
            // std::cout << "KPAIR DIFFERENT: " << first_kmer << std::endl << "KPAIR DIFFERENT: " << second_kmer << std::endl;
            return first_kmer < second_kmer;
        }

        // 4 - If equals => true if second skmer is the first one to miss a nucleotide (left based)
        else if (first_missing_nucl != second_missing_nucl){
            // std::cout << "first_missing_nucl: " << first_missing_nucl << "; second_missing_nucl: " << second_missing_nucl << std::endl;
            return first_missing_nucl < second_missing_nucl;
        }

        // 5 - If have the same position, I must compare them on all their nucleotides
        else{
            // std::cout << "same kpair, same position. " << std::endl;
            return kmer_compare(first_skmer, second_skmer, first_kmer_pos) < 0;
        }
    }

    // /** Compare 2 kmers included in 2 skmers.
    //  * @param first_skmer First kmer is included in this skmer
    //  * @param first_kmer_pos Position of the fist kmer in the first skmer. 0 is the kmer that contains the whole prefix.
    //  * @param second_skmer Second kmer is included in this skmer
    //  * @param second_kmer_pos Position of the second kmer in the second skmer. 0 is the kmer that contains the whole prefix.
    //  * @return true if the first kmer is less than the second one
    //  **/
    // bool inline kmer_less_than_kmer(const Skmer<kuint>& first_skmer, const Skmer<kuint>& second_skmer, const uint64_t kmer_pos) const
    // {
    //     auto first_kmer {first_skmer.m_pair};
    //     auto second_kmer {second_skmer.m_pair};
    //     first_kmer &= kmer_masks[kmer_pos];
    //     second_kmer &= kmer_masks[kmer_pos];
    //     return first_kmer < second_kmer;
    // }

    /** Compare 2 kmers included in 2 skmers.
     * @param first_skmer First kmer is included in this skmer
     * @param first_kmer_pos Position of the fist kmer in the first skmer. 0 is the kmer that contains the whole prefix.
     * @param second_skmer Second kmer is included in this skmer
     * @param second_kmer_pos Position of the second kmer in the second skmer. 0 is the kmer that contains the whole prefix.
     * @return true if the first kmer is less than the second one
     **/
    inline int kmer_compare(const Skmer<kuint>& first_skmer, const Skmer<kuint>& second_skmer, const uint64_t kmer_pos) const
    {
        auto first_kmer {first_skmer.m_pair};
        auto second_kmer {second_skmer.m_pair};
        first_kmer &= kmer_masks[kmer_pos];
        second_kmer &= kmer_masks[kmer_pos];

        if (first_kmer < second_kmer) return -1;
        else if (first_kmer > second_kmer) return 1;
        else return 0;
    }

    /** Masked k-mer key of `skmer` at column `kmer_pos`: the interleaved pair AND the column mask —
     * exactly the value kmer_compare compares. Exposed so the wide-store merge can cache it per cursor
     * (mask once per advance instead of twice per comparison) and compare with pair::compare3. **/
    inline kpair masked_kmer(const Skmer<kuint>& skmer, const uint64_t kmer_pos) const
    {
        return skmer.m_pair & kmer_masks[kmer_pos];
    }

    /** Check if a skmer has a kmer starting at the given position.
     * @param skmer The skmer you want to evaluate having a kmer at the given position
     * @param kmer_pos Position of the start of the kmer
     * @return true if the skmer has a valid kmer at the given position, false otherwise
     **/
    bool inline has_valid_kmer(const Skmer<kuint>& skmer, const uint64_t kmer_pos) const {
        assert(kmer_pos <= this->k - this->m);
        if ((this->m_pref_size - (this->m + 1)/2 - kmer_pos <= skmer.m_pref_size) && (skmer.m_suff_size >= kmer_pos)){
            return true;
        }
        else return false;
    }

    std::pair<uint64_t, uint64_t> get_valid_kmer_bounds(const Skmer<kuint>& skmer) const {
        const uint64_t start_pos {this->m_pref_size - skmer.m_pref_size - (this->m + 1) / 2};
        const uint64_t end_pos {skmer.m_suff_size}; //this->k - this->m - this->m_suff_size +
        return {start_pos, end_pos};

    }


    /** Generate a skmer has a kmer starting at the given position.
     * @param given_skmer The skmer from which you want to extract the kmer
     * @param kmer_pos Position of the start of the kmer
     * @return true if the skmer has a valid kmer at the given position, false otherwise
     **/
    Skmer<kuint> get_skmer_of_kmer(Skmer<kuint> given_skmer, uint64_t kmer_pos) {
        // prefix and suffix sizes computation
        // uint64_t const half_size = (2 * this->k - this->m + 1) / 2;

        // Prefix size: how many nucleotides are in the first half of the skmer
        const uint16_t prefix_size = this->m_pref_size - (this->m + 1)/2 - kmer_pos;

        // Suffix size: how many nucleotides nucleotides are in the second half of the skmer
        const uint16_t suffix_size = kmer_pos;//(this->k - prefix_size);

        // getting the kmer and generating the skmer
        kpair kmer = extract_kmer(given_skmer, kmer_pos); // extracting the kpair from the kmer

        kpair mmask = this->m_mask;

        // std::cerr << "K: " << this->k << std::endl;
        // std::cerr << "M: " << this->m << std::endl;
        // std::cerr << "GIVEN KMER: " << given_skmer << std::endl;
        // std::cerr << "[get_skmer_of_kmer] KMER POS: " << kmer_pos << std::endl;
        // std::cerr << "EXTRACTED KMER: " << kmer << std::endl;
        // std::cerr << " M_MASK: " << mmask << std::endl;
        // std::cerr << " ~KMER_MASK: " << ~kmer_masks[kmer_pos] << std::endl;
        // std::cerr << "[get_skmer_of_kmer] PREFIX SIZE: " << prefix_size << std::endl;
        // std::cerr << "[get_skmer_of_kmer] SUFFIX SIZE: " << suffix_size << std::endl;

        kmer |= (~kmer_masks[kmer_pos] & mmask); // setting to 1s the positions not used in the skmer
        // std::cerr << "OUT_SKMER: " << kmer << std::endl;
        Skmer<kuint> new_sorted_skmer(kmer, prefix_size, suffix_size);

        return new_sorted_skmer;
    }

    void inline clean_nucleotide_position_skmer(Skmer<kuint> & given_skmer, uint64_t kmer_pos) const {
        given_skmer.m_pair.m_value[0] &= (~nucleotide_masks[kmer_pos].m_value[0] & this->m_mask.m_value[0]);
        given_skmer.m_pair.m_value[1] &= (~nucleotide_masks[kmer_pos].m_value[1] & this->m_mask.m_value[1]);
        return;
    }

    std::vector<kpair > get_sp_mask()
    {
        return this->prefix_suffix_mask;
    }
    std::vector<kpair > get_k_mask()
    {
        return this->kmer_masks;
    }
    std::vector<kpair > get_n_mask()
    {
        return this->nucleotide_masks;
    }

    std::vector<kpair > generate_masks_sp()
    {
        // generate first empty SKmer
        const kuint keep_nucl = 0b11U;
        const kuint discard_nucl = 0b00U;

        for (uint64_t position {0}; position < (k - 1); position+=1){
            add_nucleotide(keep_nucl);
        }

        std::vector<kpair > masks(k - m + 2);

        for (int64_t position {static_cast<int64_t>(k - m + 1)}; position >= 0; position -= 1){
            //filling the array
            masks[position] = m_fwd.m_pair;
            //then start adding 00s while returning each time
            add_nucleotide(discard_nucl);
        }
        this->init_skmer();
        // return the vector
        return masks;
    }

    std::vector<kpair > generate_masks_k()
    {
        // generate first empty SKmer
        const kuint keep_nucl = 0b11U;
        const kuint discard_nucl = 0b00U;

        for (uint64_t position {0}; position < k; position+=1){
            add_nucleotide(keep_nucl);
        }

        std::vector<kpair > masks(k - m + 1);

        for (int64_t position {static_cast<int64_t>(k - m)}; position >= 0; position -= 1){
            //filling the array
            masks[position] = m_fwd.m_pair;
            //then start adding 00s while returning each time
            add_nucleotide(discard_nucl);
        }
        this->init_skmer();
        // return the vector
        return masks;
    }

    std::vector<kpair > generate_masks_nucleotide()
    {
        // generate first empty SKmer
        const kuint keep_nucl = 0b11U;
        const kuint discard_nucl = 0b00U;

        add_nucleotide(keep_nucl);

        std::vector<kpair > masks(2 * k - m);

        for (int64_t position {static_cast<int64_t>(2 * k - m - 1)}; position >= 0; position -= 1){
            //filling the array
            masks[position] = m_fwd.m_pair;
            //then start adding 00s while returning each time
            add_nucleotide(discard_nucl);
        }
        this->init_skmer();
        // return the vector
        return masks;
    }

    /** Rebuild the masks the query path relies on for a quotiented (b>0) manipulator. The generators
     * above run add_nucleotide with the narrowed m_mask, which wraps the suffix buffer at the
     * truncated boundary and corrupts every mask. A full-width manipulator produces the correct
     * masks; we truncate them to this store width. The reference must hold the *whole* (un-quotiented)
     * skmer, 2*(2k-m) bits. Since k is a runtime value, we cannot derive the needed width from this
     * manipulator's compile-time store type (e.g. a 258-bit skmer quotients down to a __uint128_t
     * store, but a __uint128_t pair holds only 256 bits), so we use the widest backend kuint256
     * (512-bit pair) — it holds every supported skmer and yields identical low bits for narrow stores.
     * kmer_masks drive kmer_compare; the other two are rebuilt for consistency.
     **/
    void regenerate_quotient_masks()
    {
        using full_t = kuint256;
        const uint64_t eff {2 * sk_size - m_quotient_bits};
        SkmerManipulator<full_t> full(k, m); // b == 0 -> correct full-width masks, no recursion
        const std::vector<typename Skmer<full_t>::pair> kfull  {full.get_k_mask()};
        const std::vector<typename Skmer<full_t>::pair> spfull {full.get_sp_mask()};
        const std::vector<typename Skmer<full_t>::pair> nfull  {full.get_n_mask()};
        for (size_t i {0}; i < kmer_masks.size(); i++)
            kmer_masks[i] = truncate_pair<full_t, kuint>(kfull[i], eff);
        for (size_t i {0}; i < prefix_suffix_mask.size(); i++)
            prefix_suffix_mask[i] = truncate_pair<full_t, kuint>(spfull[i], eff);
        for (size_t i {0}; i < nucleotide_masks.size(); i++)
            nucleotide_masks[i] = truncate_pair<full_t, kuint>(nfull[i], eff);
    }

    /** Returns the (k-1)-mer (prefix/suffix) starting at the given position.
     * @param skmer The skmer you want to evaluate having a kmer at the given position
     * @param start_pos Position of the start of the kmer
     * @return the k_pair associated to the k-1 mer
     **/
    kpair extract_prefix_suffix(const Skmer<kuint>& skmer, const uint64_t start_pos) {
        return skmer.m_pair & prefix_suffix_mask[start_pos];
    }

    /** Returns the (k-1)-mer (prefix/suffix) starting at the given position.
     * @param skmer The skmer you want to evaluate having a kmer at the given position
     * @param start_pos Position of the start of the kmer
     * @return the k_pair associated to the k-1 mer
     **/
    kpair extract_kmer(const Skmer<kuint>& skmer, const uint64_t start_pos) {
        assert(start_pos <= this->k - this->m);
        return skmer.m_pair & kmer_masks[start_pos];
    }

    kpair inline extract_nucleotide(const Skmer<kuint>& skmer, const uint64_t pos)
    {
        return skmer.m_pair & nucleotide_masks[pos];
    }

    /** Concatenate a kmer(given as superkmer) into an existing superkmer
    * @param skmer The skmer you want to concatenate the kmer to
    * @param kmer_skmer The skmer containing the kmer you want to concatenate
    **/
    void concatenate_skmer(Skmer<kuint>& skmer, const Skmer<kuint> kmer_skmer) const
    {
        skmer.m_pair &= kmer_skmer.m_pair;
        skmer.m_suff_size += 1;
        // sk_size(2*k-m), m_suff_size(sk_size / 2), m_pref_size((sk_size+1) / 2)
    }

    template<typename T>
    friend std::ostream& operator<<(std::ostream& os, SkmerManipulator<T>& manip);

};


template<typename T>
std::ostream& operator<<(std::ostream& os, SkmerManipulator<T>& manip)
{
    static const char nucleotides[] = {'A', 'C', 'T', 'G'};

    // cout << manip.m_fwd_suffix_buff << endl;
    // cout << manip.m_fwd_prefix_buff << endl;

    os << "[not interleaved: ";

    // Forward prefix
    for (uint64_t pref_idx{0} ; pref_idx<manip.m_pref_size ; pref_idx++)
    {
        os << (nucleotides[(manip.m_fwd_prefix_buff >> (4 * pref_idx)) & 0b11UL]);
    }
    os << " ";
    // Forward suffix
    for (uint64_t suf_idx{manip.m_suff_size} ; suf_idx>0 ; suf_idx--)
    {
        os << (nucleotides[(manip.m_fwd_suffix_buff >> (4 * suf_idx - 2)) & 0b11UL]);
    }

    os << " / ";

    // Reverse prefix
    for (uint64_t pref_idx{0} ; pref_idx<manip.m_pref_size ; pref_idx++)
    {
        os << nucleotides[(manip.m_rev_prefix_buff >> (4 * pref_idx)) & 0b11UL];
    }
    os << " ";
    // Reverse suffix
    for (uint64_t suf_idx{manip.m_suff_size} ; suf_idx>0 ; suf_idx--)
    {
        os << nucleotides[(manip.m_rev_suffix_buff >> (4 * suf_idx - 2)) & 0b11UL];
    }

    os << "]";

    return os;
}


/** Down-convert a full-width (`gen`) super-k-mer to the narrower storage type (`store`) of a
 * quotiented (b>0) sorted list. Keeps only the low `2*(2k-m) - b` bits of the interleaved pair:
 * the top `b` bits of the φ-minimizer are implied by the bucket id and dropped. Prefix/suffix
 * sizes carry over. With `store == gen` and `b == 0` this is a width-preserving copy. The caller
 * must size `store` so its pair holds the retained `2*(2k-m) - b` bits.
 *
 * Within one bucket every skmer shares the dropped top `b` bits, so truncation preserves all
 * intra-bucket k-mer comparisons and prefix/suffix extractions — query results are unchanged.
 **/
template<typename gen, typename store>
inline Skmer<store> truncate_skmer(uint64_t k, uint64_t m, uint64_t b, const Skmer<gen>& s)
{
    const uint64_t eff {2 * (2 * k - m) - b};   // retained bits
    return Skmer<store>(truncate_pair<gen, store>(s.m_pair, eff), s.m_pref_size, s.m_suff_size);
}

};

#endif