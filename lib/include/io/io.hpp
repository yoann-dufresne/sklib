#include <io/Skmer.hpp>
#include <io/Skmerator.hpp>
#include <algorithms/SkmerSorting.hpp>

#ifndef IO_H
#define IO_H


namespace km
{

template<typename kuint>
using LList = std::forward_list<Virtual_skmer<kuint>>;

template<typename kuint>
void dumpChunk(const sorting::Virtual_skmer<kuint>* dump_vector, uint64_t count, std::ofstream& outFile) {
    if (count == 0) return ;
    outFile.write(reinterpret_cast<const char*>(dump_vector), count * sizeof(sorting::Virtual_skmer<kuint>>));
    if (!outFile) {
    std::cerr << "Error writing to file. Check that the path exists." << std::endl;
    return;
    }
}


template<typename kuint>
int write_sorted_list(const SkmerManipulator<kuint> & m_manip, const LList<kuint> & list, size_t megabytes, std::string output_filename) {

    uint64_t num_bytes_per_virtual_skmer = sizeof(sorting::Virtual_skmer<kuint>>); // # of bytes per object
    uint64_t max_num_elements = megabytes / num_bytes_per_virtual_skmer;

    auto iterator = list.begin();
    int count {0};

    // 1 - Counting number of elements in the sorted list 
    while (iterator != list.end()){
        iterator++;
        count++;
    }

    // If the entire list is smaller than the chunk threshold, use the count
    uint64_t num_elements_per_dump = (count < max_num_elements) ? count : max_num_elements ;
    uint64_t num_dump_loops = count / num_elements_per_dump ;
    uint64_t remaining_elements_to_dump = count % num_elements_per_dump ;

    // 2 - setting the outstream
    std::ofstream outFile(output_filename, std::ios::binary);
    if (!outFile) { 
        std::cerr << "Error opening file. Could not store data." << std::endl;
        return 1;
    }

    // 3 - First dumping k, m and the # of elements in the vector
    outFile.write(reinterpret_cast<const char*>(&m_manip.k), sizeof(uint64_t));
    outFile.write(reinterpret_cast<const char*>(&m_manip.m), sizeof(uint64_t));
    outFile.write(reinterpret_cast<const char*>(&count), sizeof(uint64_t));
    if (!outFile) {
    std::cerr << "Error writing to file. Check that the path exists." << std::endl;
    return 1;
    }

    // 4 - Initilaizing dump vector used to store the sorted virtual skmer list
    std::vector <sorting::Virtual_skmer<kuint>> dump_vector;
    dump_vector.reserve(num_elements_per_dump);
    iterator = list.begin();
    

    // 5 - dumping 'num_elements' elements at a time, in order to not use too much memory in the output operation, as it is stored as a vector
    for (int i = 0; i < num_dump_loops; i++){
        dump_vector.clear();
        
        for (int j = 0; j < num_elements_per_dump; j++){
            dump_vector.push_back(*iterator) ;
            iterator++;
        }

        dumpChunk(dump_vector.data(), num_elements_per_dump, outFile );

        if (!outFile) {
            std::cerr << "Error occurred when writing the chunk " << i << " to file." << std::endl;
            return 1;
        }

    }

    // final dump if there are > 0 remaining elements to save
    if (remaining_elements_to_dump > 0){
        dump_vector.clear();
        dump_vector.reserve(remaining_elements_to_dump);
        for (int j = 0; j < remaining_elements_to_dump; j++){
            dump_vector.push_back(*iterator) ;
            iterator++;
        }
        dumpChunk(dump_vector.data(), remaining_elements_to_dump, outFile );
        if (!outFile) {
            std::cerr << "Error occurred when writing the chunk " << i << " to file." << std::endl;
            return 1;
        }

    }

    outFile.close();
    return 0;
}

template<typename kuint>
int load_sorted_vector(km::SkmerManipulator<kuint> & m_manip, std::vector <sorting::Virtual_skmer<kuint>> & sorted_vector, std::string input_filename){

    // 1 - setting the input stream
    std::ifstream inFile(input_filename, std::ios::binary);

    if (!input_filename) { 
        std::cerr << "Error opening input file. Could not read." << std::endl;
        return 1;
    }

    // 2 - reading the first 3 parameters
    uint64_t count, k, m ;
    inFile.read(reinterpret_cast<char*>(&k), sizeof(uint64_t));
    inFile.read(reinterpret_cast<char*>(&m), sizeof(uint64_t));
    inFile.read(reinterpret_cast<char*>(&count), sizeof(uint64_t));
    if (!inFile) {
        std::cerr << "Error reading header from file." << std::endl;
    return 1;
    }

    // 3 - load into memory the Virtual Super Kmer vector
    sorted_vector.resize(count);
    inFile.read(reinterpret_cast<char*>(sorted_vector.data()), count * sizeof(sorting::Virtual_skmer<kuint>));
    if (!inFile) {
        std::cerr << "Error reading Virtual_skmer data from file." << std::endl;
        return 1;
    }
    inFile.close();

    // 4 - generate the manipulator
    &m_manip = SkmerManipulator(k, m);
    
    return 0;

}

}