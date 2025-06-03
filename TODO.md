### FINISHING PROJECT TODO
- [X] WRITE SORTED SKMER VECTOR DISK DUMP
- [X] WRITE SORTED SKMER VECTOR DISK LOAD
- [X] TEST DISK DUMP / LOAD
- [ ] WIRTE CLASS WRAPPER OF SORTED SKMER LIST
- [ ] USE THE WRAPPER CLASS FOR DISK LOAD AND DUMP
- [ ] QUERY A SINGLE KMER FROM THE LIST
- [ ] QUERY A STREAM OF KMERS FROM THE LIST


#### OLD TODOs, TO CHECK
* Add the Bzip2 dependancy in the CMake
* In Skmer header, Check how "pair& operator|= (const kuint& value)" behaves related to "pair(kuint& single) : m_value(0, single)"
* when adding a single kuint, should go to least significant not most.