#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/SkmerSorting.hpp>

#include <fstream>

#ifndef IO_H
#define IO_H

constexpr uint64_t num_bytes_in_megabyte = 1000000 ;

namespace km
{
namespace io
{
template<typename kuint>
using LList = std::forward_list<km::sorting::Virtual_skmer<kuint>>;

template<typename kuint>
void dumpChunk(const Skmer<kuint>* dump_vector, uint64_t count, std::ofstream& outFile) {
    if (count == 0) return ;
    outFile.write(reinterpret_cast<const char*>(dump_vector), count * sizeof(Skmer<kuint>));
    if ( outFile.fail() ) {
    std::cerr << "Error writing to file. Check that the path exists." << std::endl;
    return;
    }
}


template<typename kuint>
int write_sorted_list(const SkmerManipulator<kuint> & m_manip, const LList<kuint> & list, size_t megabytes, std::string output_filename) {

    uint64_t num_bytes_per_virtual_skmer = sizeof(sorting::Virtual_skmer<kuint>); // # of bytes per object
    
    // std::cerr << "MEGABYTES: " << megabytes << " , NUM_BYTES_PER_VSKMER: " << num_bytes_per_virtual_skmer << " ." << std::endl;
    
    uint64_t max_num_elements = megabytes * num_bytes_in_megabyte / num_bytes_per_virtual_skmer;
    
    // std::cerr << "MAX NUM ELEMENTS: " << max_num_elements << std::endl;

    auto iterator = list.begin();
    uint64_t count {0};

    // 1 - Counting number of elements in the sorted list 
    while (iterator != list.end()){
        iterator++;
        count++;
    }
    // std::cerr << "NUM ELEMENTS IN LIST: " << count << std::endl;

    // If the entire list is smaller than the chunk threshold, use the count
    uint64_t num_elements_per_dump = (count < max_num_elements) ? count : max_num_elements ;
    uint64_t num_dump_loops = count / num_elements_per_dump ;
    uint64_t remaining_elements_to_dump = count % num_elements_per_dump ;

    // std::cerr << "NUM ELEMENTS PER DUMP: " << num_elements_per_dump << std::endl;
    // std::cerr << "NUM DUMP LOOPS: " << num_dump_loops << std::endl;
    // std::cerr << "NUM ELEMENTS LAST DUMP: " << remaining_elements_to_dump << std::endl;

    // 2 - setting the outstream
    std::ofstream outFile(output_filename, std::ios::binary);
    if ( outFile.fail()) { 
        std::cerr << "Error opening file. Could not store data." << std::endl;
        return 1;
    }

    // std::cerr << "WRITING TO DISK..." << std::endl;
    // std::cerr << "k: " << m_manip.k << std::endl;
    // std::cerr << "m: " << m_manip.m << std::endl;
    // std::cerr << "count: " << count << std::endl;

    // 3 - First dumping k, m and the # of elements in the vector
    outFile.write(reinterpret_cast<const char*>(&m_manip.k), sizeof(uint64_t));
    outFile.write(reinterpret_cast<const char*>(&m_manip.m), sizeof(uint64_t));
    outFile.write(reinterpret_cast<const char*>(&count), sizeof(uint64_t));
    if (!outFile) {
    std::cerr << "Error writing to file. Check that the path exists." << std::endl;
    return 1;
    }

    // std::cerr << "FINISHED WRITING VARIABLES TO DISK." << std::endl;

    // 4 - Initilaizing dump vector used to store the sorted virtual skmer list
    std::vector <Skmer<kuint>> dump_vector;
    dump_vector.reserve(num_elements_per_dump);
    iterator = list.begin();
    

    // 5 - dumping 'num_elements' elements at a time, in order to not use too much memory in the output operation, as it is stored as a vector
    for (int i = 0; i < num_dump_loops; i++){
        dump_vector.clear();
        
        for (int j = 0; j < num_elements_per_dump; j++){
            dump_vector.push_back((*iterator).skmer) ;
            iterator++;
        }

        // std::cerr << "DUMPING A CHUNK OF DATA." << std::endl;
        dumpChunk(dump_vector.data(), num_elements_per_dump, outFile );

        if (!outFile ) {
            std::cerr << "Error occurred when writing the chunk " << i << " to file." << std::endl;
            return 1;
        }

    }

    // final dump if there are > 0 remaining elements to save
    if (remaining_elements_to_dump > 0){
        dump_vector.clear();
        dump_vector.reserve(remaining_elements_to_dump);
        for (int j = 0; j < remaining_elements_to_dump; j++){
            dump_vector.push_back((*iterator).skmer) ;
            iterator++;
        }
        dumpChunk(dump_vector.data(), remaining_elements_to_dump, outFile );
        if (!outFile) {
            std::cerr << "Error occurred when writing the last chunk to file." << std::endl;
            return 1;
        }

    }

    // std::cerr << "FINISHED WRITING TO DISK." << std::endl;
    outFile.close();
    return 0;
}

template<typename kuint>
void load_sorted_vector(const std::string input_filename, uint64_t& m_k, uint64_t& m_m, std::vector<Skmer<kuint>>& sorted_vector){

    // 1 - setting the input stream
    // std::cerr << "SETTING INPUT STREAM." << std::endl;
    std::ifstream inFile(input_filename, std::ios::binary);

    if (inFile.fail()) { 
        std::cerr << "Error opening input file. Could not read." << std::endl;
        throw std::runtime_error("Failed to load sorted vector from " + input_filename);
    }

    // 2 - reading the first 3 parameters
    uint64_t count;
    inFile.read(reinterpret_cast<char*>(&m_k), sizeof(uint64_t));
    inFile.read(reinterpret_cast<char*>(&m_m), sizeof(uint64_t));
    inFile.read(reinterpret_cast<char*>(&count), sizeof(uint64_t));
    if (inFile.fail()) {
        std::cerr << "Error reading header from file." << std::endl;
        throw std::runtime_error("Failed to load sorted vector from " + input_filename);
    }

    // std::cerr << "FINISHED READING VARIABLES FROM DISK." << std::endl;
    // std::cerr << "k: " << m_k << std::endl;
    // std::cerr << "M: " << m_m << std::endl;
    // std::cerr << "count: " << count << std::endl;

    // 3 - load into memory the Virtual Super Kmer vector
    // std::vector<Skmer<kuint>> sorted_vector;
    sorted_vector.reserve(count);
    inFile.read(reinterpret_cast<char*>(sorted_vector.data()), count * sizeof(Skmer<kuint>));
    if (inFile.fail()) {
        std::cerr << "Error reading Virtual_skmer data from file." << std::endl;
        throw std::runtime_error("Failed to load sorted vector from " + input_filename);
    }
    inFile.close();

    // std::cerr << "FINISHED READING FROM DISK." << std::endl;
    return;

}

}
}

#endif