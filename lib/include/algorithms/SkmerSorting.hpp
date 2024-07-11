#include <array>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <unordered_map>

#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>

#ifndef SKMERSORTING_H
#define SKMERSORTING_H


namespace km
{
    namespace sorting
    {

    template <class It, typename kuint>
    class compare_kmer_skmer_pos;


    /** Sorts skmer ids based on the kmers they contain at a given positon.
     * @param start start_position in the skmer generator
     * @param end end_positon in the skmer generator
     * @param kmer_pos position of the kmer in the skmer (column position)
     * @param m_manip skmer manipulator
     * @return a vector of sorted superkmer ids (if no kmer, no skmer id) 
     **/
    template<class It, typename kuint>
    std::vector<uint64_t> sort_column(It start, It end, uint64_t kmer_pos, SkmerManipulator<kuint>& m_manip)
    {   
        // Accessing and comparing kmers in skmers (less than) is done by kmer_lt_kmer of skmermanipulator
        // 1st pass over the column: check which skmers are ok to be processed
        // Check if the first skmer has a kmer in this position
        std::vector<uint64_t> valid_skmer;
        uint64_t sk_id = 0;
        km::SkmerPrettyPrinter<kuint> pp {m_manip.k, m_manip.m};
        //Iterating over the range [start, end)
        for(It skmer = start; skmer != end; ++skmer)
        {
            pp << *skmer;
            std::cout << "checking kmer validity" << pp << std::endl;
            if (m_manip.has_valid_kmer(*skmer, kmer_pos)){
                valid_skmer.push_back(sk_id);
                std::cout << "valid" << std::endl;
            }
            sk_id++;
        }

        // 2nd pass over the column: return ordered list 
        // For every "column" i.e. possible kmer in the skmer size
        // For every skmer that has a kmer in that column
        std::sort(valid_skmer.begin(), valid_skmer.end(),
                compare_kmer_skmer_pos<It, kuint>(kmer_pos, m_manip, start, end));

        std::cout << "SORTED SKMER LIST - ( size: " << valid_skmer.size() << ") " << std::endl;
        for (uint64_t i: valid_skmer) 
            std::cout << i << ' ';
        std::cout << std::endl;
        
        return valid_skmer;
    }

    /** Returns candidate overlaps between two columns of sorted skmer ids
     * @param left_column first column
     * @param right_column second column (contigous)
     * @return a vector of pairs of candidate overlaps between the two columns
     **/
    template<class It, typename kuint>
    std::vector<std::pair<uint64_t, uint64_t> > get_candidate_overlaps(std::vector<Skmer<kuint> > & skmer_enumeration, SkmerManipulator<kuint>& manipulator, uint64_t left_position, std::vector<uint64_t> left_column, std::vector<uint64_t> right_column)
    {
        std::unordered_map< template Skmer<kuint>::pair, std::vector<uint64_t> > prefixes;

        Skmer<kuint>::pair suffix, prefix;
        std::vector<std::pair<uint64_t,uint64_t> > candidare_overlaps;
        std::unordered_map<Skmer<kuint>::pair, std::vector<uint64_t> >::const_iterator matching_prefix;
        // First, there should be a function that extracts the k-1 suffix of the left column
        for (auto& skmer_id : right_column) {
            prefix = manip.extract_fix(skmer_enumeration[skmer_id], left_position);
            prefixes[prefix].push_back(skmer_id);
        }

        // Second, there should be a function that extracts the k-1 prefix of the right column (same funct as before, just give param the place)
        for (auto& skmer_id : left_column) {
            suffix = manip.extract_fix(skmer_enumeration[skmer_id], left_position+1);
            matching_prefix = prefixes.find (suffix);
            if (matching_prefix != prefixes.end()){
                for (auto& pref_sk_id: matching_prefix->second){
                    candidare_overlaps.emplace_back(skmer_id,pref_sk_id);
                }
            }
        }
        return candidare_overlaps;
    }


    template <class It, typename kuint>
    class compare_kmer_skmer_pos {
        uint64_t position;
        SkmerManipulator<kuint> & manipulator;
        const It start;
        const It end;


    public:
        // the comparison function takes as argument 2 integers, a position and the vector of skmers. 
        // It compares the two skmers in the selected position and returns which one is before the other.
        compare_kmer_skmer_pos(uint64_t p, SkmerManipulator<kuint> & skmer_manipulator, const It start_skmer_en, const It end_skmer_en) // 
        : position(p), manipulator(skmer_manipulator), start(start_skmer_en), end(end_skmer_en) {}

        bool operator()(const uint64_t skmer_id_1,const uint64_t skmer_id_2) const {
            assert((start+skmer_id_1) < end);
            assert((start+skmer_id_2) < end);
            return manipulator.kmer_lt_kmer(*(start+skmer_id_1), position, *(start+skmer_id_2), position);
        }
    };
    }
}

#endif