### FINISHING PROJECT TODO
- [X] WRITE SORTED SKMER VECTOR DISK DUMP
- [X] WRITE SORTED SKMER VECTOR DISK LOAD
- [X] TEST DISK DUMP / LOAD
- [X] BIG/LITLLE ENDIAN SANITY CHECK
- [X] WIRTE CLASS WRAPPER OF SORTED SKMER LIST
- [X] USE THE WRAPPER CLASS FOR DISK LOAD AND DUMP
- [X] REIMPLEMENT TESTS FOR THE NEW VIRTUALSKMER WRAPPER CLASS
- [X] QUERY A SINGLE KMER FROM THE LIST
- [X] QUERY A STREAM OF KMERS FROM THE LIST
- [X] UPDATE QUERY STRATEGY FOR CACHE LOCALITY (bucketed routing + lazy per-bucket load)
- [ ] ADD UNIQUENESS AFTER K-MER COLUMN SORTING
- [ ] PARELLALIZE QUERY ON MULTIPLE CORES (deferred — sklib is sequential by design for now)
- [ ] CI/CD WITH GITHUB ACTIONS
- [ ] DIFF BETWEEN 2 PRECOMPUTED LISTS
- [ ] MERGE OF 2 PRECOMPUTED LISTS
- [ ] INTERSECTION OF 2 PRECOMPUTED LISTS
- [ ] ADD HASHING
- [ ] TEST ON HUMAN, METAGENOME, BACTERIA, RICE PANGENOME?
- [X] TEST AGAINST BQF, SSHASH (benchmark vs bqf/sshash/sbwt/cbl — see scripts/bench)


#### OLD TODOs, TO CHECK
* Add the Bzip2 dependancy in the CMake
* In Skmer header, Check how "pair& operator|= (const kuint& value)" behaves related to "pair(kuint& single) : m_value(0, single)"
* when adding a single kuint, should go to least significant not most.