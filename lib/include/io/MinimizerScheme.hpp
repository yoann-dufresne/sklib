#include <cstdint>
#include <algorithm>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef SKLIB_MINIMIZER_SCHEME_HPP
#define SKLIB_MINIMIZER_SCHEME_HPP

// Minimizer sampling scheme: WHICH m-mer of a k-mer window is the minimizer.
//
// Two schemes, both canonical (strand-invariant) and both "argmin of a per-m-mer rank":
//
//   * random      — rank(x) = φ(canonical interleaved code of x). The historical scheme; the default.
//   * open-closed — rank(x) = (class(x) << 2m) | φ(canonical interleaved code of x), where class(x)
//                   is the canonical open-closed syncmer class of the m-mer (0 = open, 1 = closed,
//                   2 = other). For random, class ≡ 0, so the two ranks agree bit for bit.
//
// The class is a SELECTION key only: it is never stored. The on-disk slot stays ψ(min), so the
// record layout, the per-column order, the bucket routing and the set-op merge are untouched; only
// which records exist changes. The rank stays injective on canonical m-mers (φ is), so "equal rank
// <=> same canonical m-mer" — the premise of minimizer_is_ambiguous, the most-central tie-break,
// choose_kmer_minimizer and canonical_pieces — still holds.
//
// Open-closed class of an m-mer x (linear nucleotides x[0..m)), with o = (m - s) / 2:
//
//     h_s(u)  = π_s(min(u, rc(u)))                     canonical s-mer key
//     M       = min_{i in [0, m-s]} h_s(x[i .. i+s))
//     class   = 0  if h_s(x[o ..]) == M                 (open: the CENTRAL s-mer attains the min)
//             = 1  elif h_s(x[0 ..]) == M or h_s(x[m-s ..]) == M   (closed: an END s-mer does)
//             = 2  otherwise
//
// "Attains the minimum", not "leftmost argmin": the test must be a symmetric function of the key
// multiset. Under reverse complement the s-mer at offset i moves to offset m-s-i with the same
// canonical key, the ends {0, m-s} map to themselves, and o maps to itself iff m - s is even. Hence
// s ≡ m (mod 2) is MANDATORY: with the wrong parity the open class is not strand-symmetric and a
// k-mer stored from one strand is framed differently when queried from the other (false negative).
// See notes/minimizers/04-canonical-strand-invariant.md §6.3b.
//
// Two implementations of the class share the same key function h_s:
//   * OcClassifier::class_of — the standalone O(m) reference (ambiguous path, tests);
//   * OcStream               — the O(1)-amortized streaming version fed one base at a time by the
//                              fast path (SeqSkmerator). It is tested equal to class_of at every
//                              position (tests/km/minimizer_scheme.cpp).

namespace km
{

enum class MinimizerScheme : uint8_t { random = 0, open_closed = 1 };

// Scheme + its parameter. `s` is meaningful for open-closed only (0 for random). Recorded in the
// VSKMER_7 header as scheme_id(8) + scheme_param(8).
struct SchemeParams
{
    MinimizerScheme scheme {MinimizerScheme::random};
    uint8_t s {0};

    // s-mer keys are 64-bit words (2s bits): s <= 32. The automatic default is <= 5.
    static constexpr uint64_t MAX_OC_S {32};

    bool is_open_closed() const { return scheme == MinimizerScheme::open_closed; }

    // Bits the class adds above the 2m-bit φ value in the selection rank (the generation word must
    // hold 2m + class_bits()).
    uint64_t class_bits() const { return is_open_closed() ? 2 : 0; }

    // Automatic s for open-closed: 3 for odd m <= 17, 5 for odd m >= 19, 4 for even m >= 8, 2 for
    // even m < 8 (the even-m values are unmeasured). Clamped down, keeping the parity, to the
    // largest s <= m - 2 (so m = 3 gets s = 1). Returns 0 when no valid s exists (m <= 2).
    static uint64_t default_s(uint64_t m)
    {
        uint64_t s {(m % 2 == 1) ? (m <= 17 ? 3u : 5u) : (m >= 8 ? 4u : 2u)};
        while (s >= 1 && s + 2 > m) {
            if (s <= 2) return 0;
            s -= 2;
        }
        return s;
    }

    static SchemeParams random_scheme() { return SchemeParams{}; }

    // Open-closed with the given s, or the automatic default when s == 0. Not validated.
    static SchemeParams open_closed(uint64_t m, uint64_t s = 0)
    {
        SchemeParams p;
        p.scheme = MinimizerScheme::open_closed;
        p.s = static_cast<uint8_t>(std::min<uint64_t>(s == 0 ? default_s(m) : s, 255));
        return p;
    }

    // Empty string when the parameters are valid for this m, else a human-readable reason.
    std::string validation_error(uint64_t m) const
    {
        if (!is_open_closed()) return {};
        if (m < 3)
            return "open-closed minimizers need m >= 3 (there must be an s-mer strictly inside the "
                   "m-mer, 1 <= s <= m-2); got m=" + std::to_string(m);
        if (s < 1 || s + 2 > m)
            return "--oc-s must satisfy 1 <= s <= m-2 = " + std::to_string(m - 2) +
                   " (got s=" + std::to_string(s) + ")";
        if (s > MAX_OC_S)
            return "--oc-s must be <= " + std::to_string(MAX_OC_S) + " (got s=" + std::to_string(s) + ")";
        if ((m - s) % 2 != 0) {
            const uint64_t dflt {default_s(m)};
            return "--oc-s must have the same parity as m (s ≡ m mod 2; got m=" + std::to_string(m) +
                   ", s=" + std::to_string(s) + "). The open-syncmer test looks at the s-mer at offset "
                   "(m-s)/2, which is the exact centre of the m-mer only when m-s is even; with the "
                   "other parity the class is not reverse-complement symmetric and k-mers queried from "
                   "the opposite strand would be missed." +
                   (dflt ? " For m=" + std::to_string(m) + " the default is s=" + std::to_string(dflt) + "."
                         : std::string{});
        }
        return {};
    }

    void validate(uint64_t m) const
    {
        const std::string e {validation_error(m)};
        if (!e.empty()) throw std::invalid_argument(e);
    }

    static const char* name(MinimizerScheme sch)
    {
        return sch == MinimizerScheme::open_closed ? "open-closed" : "random";
    }

    // "random" / "open-closed s=3": used in error messages and logs.
    std::string to_string() const
    {
        if (!is_open_closed()) return name(scheme);
        return std::string{name(scheme)} + " s=" + std::to_string(s);
    }

    // CLI spelling -> scheme ("random" | "open-closed").
    static std::optional<MinimizerScheme> parse(const std::string& v)
    {
        if (v == "random") return MinimizerScheme::random;
        if (v == "open-closed") return MinimizerScheme::open_closed;
        return std::nullopt;
    }

    // On-disk encoding (VSKMER_7 header).
    uint64_t header_id() const { return static_cast<uint64_t>(scheme); }
    uint64_t header_param() const { return is_open_closed() ? s : 0; }

    // Decode a VSKMER_7 header's scheme fields, rejecting unknown ids and invalid parameters (a file
    // written by a newer sklib, or corrupted).
    static SchemeParams from_header(uint64_t id, uint64_t param, uint64_t m, const std::string& filename)
    {
        SchemeParams p;
        if (id == static_cast<uint64_t>(MinimizerScheme::random)) {
            if (param != 0)
                throw std::runtime_error("Invalid minimizer-scheme parameter " + std::to_string(param) +
                                         " for the random scheme in " + filename);
            return p;
        }
        if (id == static_cast<uint64_t>(MinimizerScheme::open_closed)) {
            if (param > 255)
                throw std::runtime_error("Invalid open-closed parameter s=" + std::to_string(param) + " in " + filename);
            p.scheme = MinimizerScheme::open_closed;
            p.s = static_cast<uint8_t>(param);
            const std::string e {p.validation_error(m)};
            if (!e.empty())
                throw std::runtime_error("Invalid open-closed parameters in " + filename + ": " + e);
            return p;
        }
        throw std::runtime_error("Unknown minimizer scheme id " + std::to_string(id) + " in " + filename +
                                 " (written by a newer sklib?)");
    }

    friend bool operator==(const SchemeParams& a, const SchemeParams& b)
    {
        return a.scheme == b.scheme && a.header_param() == b.header_param();
    }
    friend bool operator!=(const SchemeParams& a, const SchemeParams& b) { return !(a == b); }
};


// Canonical s-mer keys and the reference open-closed class of an m-mer.
//
// s-mer codes are "linear": first nucleotide in the most significant bits (2-bit encoding of
// Skmer.hpp: A=0, C=1, T=2, G=3; complement = x ^ 2). The canonical code is min(code, rc code) and
// the key is π_s of it, π_s a fixed xorshift-multiply bijection on 2s bits (so keys tie only for
// occurrences of the same canonical s-mer). For s <= TABLE_MAX_S the key of every forward code is
// precomputed (4^s uint16 entries, <= 128 KiB at s = 8), which makes the streaming per-base cost a
// single table lookup.
class OcClassifier
{
public:
    static constexpr uint64_t TABLE_MAX_S {8};

    OcClassifier(uint64_t m, uint64_t s)
        : m_m(m), m_s(s), m_o((m - s) / 2)
        , m_mask(s >= 32 ? ~uint64_t{0} : ((uint64_t{1} << (2 * s)) - 1))
        , m_sh1(std::max<uint64_t>(1, s)), m_sh2(std::max<uint64_t>(1, s > 1 ? s - 1 : 1))
    {
        SchemeParams p;
        p.scheme = MinimizerScheme::open_closed;
        p.s = static_cast<uint8_t>(std::min<uint64_t>(s, 255));
        if (s != p.s || !p.validation_error(m).empty())
            throw std::invalid_argument("OcClassifier: invalid open-closed parameters (m=" + std::to_string(m) +
                                        ", s=" + std::to_string(s) + "): " + p.validation_error(m));
        if (s <= TABLE_MAX_S) {
            m_table.resize(uint64_t{1} << (2 * s));
            for (uint64_t u {0}; u < m_table.size(); ++u)
                m_table[u] = static_cast<uint16_t>(key_from_codes(u, rc_code(u, s)));
        }
    }

    uint64_t m() const { return m_m; }
    uint64_t s() const { return m_s; }
    uint64_t open_offset() const { return m_o; }               // o = (m - s) / 2
    uint64_t smers_per_mmer() const { return m_m - m_s + 1; }  // W
    uint64_t smer_mask() const { return m_mask; }
    bool has_table() const { return !m_table.empty(); }
    // forward code -> h_s, or nullptr when s > TABLE_MAX_S (then use key_from_codes).
    const uint16_t* table() const { return m_table.empty() ? nullptr : m_table.data(); }

    // Reverse complement of a linear code of `len` nucleotides.
    static uint64_t rc_code(uint64_t code, uint64_t len)
    {
        uint64_t r {0};
        for (uint64_t i {0}; i < len; ++i) {
            r = (r << 2) | ((code & 0b11) ^ 0b10);
            code >>= 2;
        }
        return r;
    }

    // π_s: bijection on 2s bits (each step is invertible: xorshift by >= 1, odd multiply mod 2^2s,
    // xor with a constant). Constants differ from φ's so the s-mer order is unrelated to the m-mer
    // order.
    uint64_t pi(uint64_t x) const
    {
        x &= m_mask;
        x ^= x >> m_sh1;
        x = (x * 0xff51afd7ed558ccdULL) & m_mask;
        x ^= x >> m_sh2;
        x = (x * 0xc4ceb9fe1a85ec53ULL) & m_mask;
        x ^= x >> m_sh1;
        x ^= 0x2545F4914F6CDD1DULL & m_mask;
        return x & m_mask;
    }

    // h_s from the forward and reverse-complement codes of the same s-mer.
    uint64_t key_from_codes(uint64_t fwd, uint64_t rc) const { return pi(std::min(fwd, rc)); }

    // h_s from both codes, through the table when there is one (the codes must be consistent).
    uint64_t key_fast(uint64_t fwd, uint64_t rc) const
    {
        return has_table() ? m_table[fwd] : key_from_codes(fwd, rc);
    }

    // h_s of the s nucleotides x[0..s).
    template<typename T>
    uint64_t key_at(const T* x) const
    {
        uint64_t f {0};
        for (uint64_t j {0}; j < m_s; ++j) f = (f << 2) | (static_cast<uint64_t>(x[j]) & 0b11);
        return key_from_codes(f, rc_code(f, m_s));
    }

    // Reference open-closed class of the m-mer x[0..m) (see the header comment). O(m·s); used by
    // the ambiguous path (memoized per position) and as the ground truth for OcStream.
    template<typename T>
    uint8_t class_of(const T* x) const
    {
        const uint64_t W {smers_per_mmer()};
        uint64_t M {~uint64_t{0}};
        for (uint64_t i {0}; i < W; ++i) M = std::min(M, key_at(x + i));
        if (key_at(x + m_o) == M) return 0;
        if (key_at(x) == M || key_at(x + (W - 1)) == M) return 1;
        return 2;
    }

private:
    uint64_t m_m, m_s, m_o, m_mask, m_sh1, m_sh2;
    std::vector<uint16_t> m_table;  // forward code -> h_s, when s <= TABLE_MAX_S
};


// Streaming open-closed class: fed every real nucleotide of a sequence in order (push), it records
// the class of each complete m-mer in a ring indexed by the m-mer's start position (class_at).
//
// Per base: roll the s-mer codes, one key (table lookup when s <= 8), and a sliding-window minimum
// over the last W = m-s+1 keys. The minimum is maintained re-scan style, like the minimizer window:
// keep the minimum value and the LAST position it occurs at; a new key <= it takes over; when that
// position leaves the window, rescan the W keys. O(1) amortized on random input, and an equal key
// refreshes the position, so homopolymers / tandem repeats never trigger rescans.
//
// The ring of classes holds `class_capacity` (>= k-m+1) consecutive m-mer starts: the Skmerator reads
// the class of the candidate's central m-mer, which ended k-m bases before the newest base.
class OcStream
{
public:
    OcStream() = default;

    void init(const OcClassifier* c, uint64_t class_capacity)
    {
        m_c = c;
        m_table = c->table();
        m_s = c->s();
        m_m = c->m();
        m_o = c->open_offset();
        m_W = c->smers_per_mmer();
        m_mask = c->smer_mask();
        m_rc_shift = 2 * (m_s - 1);
        m_keys.assign(bit_ceil_u64(m_W), 0);
        m_keys_mask = m_keys.size() - 1;
        m_cls.assign(bit_ceil_u64(std::max<uint64_t>(class_capacity, 1)), 0);
        m_cls_mask = m_cls.size() - 1;
        reset();
    }

    bool active() const { return m_c != nullptr; }

    // Start a new sequence (absolute positions restart at 0).
    void reset()
    {
        m_pos = 0;
        m_fwd = 0;
        m_rc = 0;
        m_min = ~uint64_t{0};
        m_min_pos = 0;
    }

    // Feed the next nucleotide (2-bit value). When it completes an m-mer starting at j, the class of
    // that m-mer becomes available as class_at(j).
    // (Works on locals and writes the members back once: the key ring holds uint64_t, so a store
    // into it could alias the uint64_t members and would otherwise force them to be reloaded.)
    [[gnu::always_inline]] inline void push(uint64_t nucl)
    {
        const uint64_t i {m_pos++};
        nucl &= 0b11;
        const uint64_t fwd {((m_fwd << 2) | nucl) & m_mask};
        const uint64_t rc {(m_rc >> 2) | ((nucl ^ 0b10) << m_rc_shift)};
        m_fwd = fwd;
        m_rc = rc;
        if (i + 1 < m_s) return;
        const uint64_t p {i + 1 - m_s};                      // start of the newest s-mer
        const uint64_t key {m_table ? uint64_t{m_table[fwd]} : m_c->key_from_codes(fwd, rc)};
        uint64_t* const keys {m_keys.data()};
        const uint64_t kmask {m_keys_mask};
        keys[p & kmask] = key;

        const uint64_t W {m_W};
        const uint64_t lo {(p + 1 >= W) ? p + 1 - W : 0};    // window of s-mer starts [lo, p]
        uint64_t mn {m_min}, mpos {m_min_pos};
        if (mpos < lo) {
            // The current minimum left the window: rescan, keeping the last position of the minimum.
            mn = keys[lo & kmask];
            mpos = lo;
            for (uint64_t q {lo + 1}; q <= p; ++q) {
                const uint64_t v {keys[q & kmask]};
                if (v <= mn) { mn = v; mpos = q; }
            }
        } else if (key <= mn) {
            mn = key;
            mpos = p;
        }
        m_min = mn;
        m_min_pos = mpos;

        if (p + 1 < W) return;
        const uint64_t j {lo};                                // m-mer start: its s-mers are [j, p]
        const uint8_t cls {static_cast<uint8_t>(
            (keys[(j + m_o) & kmask] == mn) ? 0 : ((keys[j & kmask] == mn || key == mn) ? 1 : 2))};
        m_cls.data()[j & m_cls_mask] = cls;
    }

    // Class of the m-mer starting at absolute position `start` (must be complete and still within the
    // ring capacity).
    uint8_t class_at(uint64_t start) const { return m_cls[start & m_cls_mask]; }

    // Number of nucleotides pushed since the last reset.
    uint64_t pushed() const { return m_pos; }

private:
    static uint64_t bit_ceil_u64(uint64_t v)
    {
        uint64_t c {1};
        while (c < v) c <<= 1;
        return c;
    }

    const OcClassifier* m_c {nullptr};
    const uint16_t* m_table {nullptr};   // the classifier's key table (s <= 8), cached
    uint64_t m_s {0}, m_m {0}, m_o {0}, m_W {0}, m_mask {0}, m_rc_shift {0};
    uint64_t m_pos {0}, m_fwd {0}, m_rc {0};
    uint64_t m_min {~uint64_t{0}}, m_min_pos {0};
    std::vector<uint64_t> m_keys;
    uint64_t m_keys_mask {0};
    std::vector<uint8_t> m_cls;
    uint64_t m_cls_mask {0};
};

} // namespace km

#endif // SKLIB_MINIMIZER_SCHEME_HPP
