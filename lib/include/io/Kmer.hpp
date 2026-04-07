#include <iostream>

#ifndef KMER_H
#define KMER_H


namespace km
{

template<typename kuint>
struct Kmer
{
	kuint m_value;

	Kmer() : Kmer(0) {};

	Kmer(const kuint value) : m_value(value)
	{}

	Kmer(const Kmer<kuint>& other) : Kmer(other.m_value)
	{}

	Kmer<kuint>& operator= (Kmer<kuint>& other)
	{ m_value = other.m_value; return *this; }
	Kmer<kuint>& operator= (Kmer<kuint>&& other)
	{ m_value = other.m_value; return *this; }

	bool operator<(const Kmer<kuint>& other) const
	{ return m_value < other.m_value; }

	bool operator==(const Kmer<kuint>& other) const
	{ return m_value == other.m_value; }

	bool operator>(const Kmer<kuint>& other) const
	{ return m_value > other.m_value; }

	template<typename T>
	friend std::ostream& operator<<(std::ostream& os, Kmer<T>& kmer);
};

template<typename T>
std::ostream& operator<<(std::ostream& os, Kmer<T>& kmer)
{
	os << kmer.m_value;
    return os;
};


		using namespace std;
template<typename kuint>
class KmerManipulator
{
public:
	const uint64_t k;
	Kmer<kuint> m_fwd;
	Kmer<kuint> m_rev;
protected:
	const kuint m_mask;
	const uint64_t m_leftmost_shift;
public:
	KmerManipulator(const uint64_t k) : k(k), m_fwd(0), m_rev(0)
		, m_mask( static_cast<kuint>(~static_cast<kuint>(0)) >> (sizeof(kuint) * 8 - 2 * k) ), m_leftmost_shift( (k-1) * 2 )
	{
		assert((k+3) / 4 <= sizeof(kuint));
	}

	void init_kmer()
	{
		m_fwd = m_rev = Kmer<kuint>{0};
	}

	void set_kmer(const Kmer<kuint>& fwd, const Kmer<kuint>& rev)
	{
		m_fwd = fwd;
		m_rev = rev;
	}

	void set_kmer(const Kmer<kuint>& fwd)
	{
		m_fwd = fwd;
		m_rev = ~static_cast<kuint>(0);
	}

	inline Kmer<kuint> add_nucleotide(const kuint nucl)
	{
		// Forward
		m_fwd.m_value >>= 2;
		m_fwd.m_value += nucl << m_leftmost_shift;

		// Reverse
		const kuint rev_nucl {static_cast<kuint>((nucl + 2U) & 0b11U)};
		m_rev.m_value <<= 2;
		m_rev.m_value &= m_mask;
		m_rev.m_value |= rev_nucl;

		return std::min(m_fwd, m_rev);
	}

	template<typename T>
	friend std::ostream& operator<<(std::ostream& os, KmerManipulator<T>& manip);
};


template<typename T>
std::ostream& operator<<(std::ostream& os, KmerManipulator<T>& manip)
{
	static const char nucleotides[] = {'A', 'C', 'T', 'G'};
	os << '[';

	// Forward
	for (uint64_t idx{0} ; idx<manip.k ; idx++)
		os << nucleotides[(manip.m_fwd.m_value >> (2 * idx)) & 0b11U];
	os << '(' << manip.m_fwd.m_value << ") / ";

	// Reverse
	for (uint64_t idx{0} ; idx<manip.k ; idx++)
		os << nucleotides[(manip.m_rev.m_value >> (2 * idx)) & 0b11U];
	os << '('  << manip.m_rev.m_value << ")]";

    return os;
}

};



#endif
