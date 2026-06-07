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
- [X] PARALLELIZE QUERY ON MULTIPLE CORES (file -i path: std::thread producer/consumer, -t/--threads default 8, output in input order, byte-identical to sequential — ParallelQuery.hpp)
- [X] BUCKETED VSKMER_4 FORMAT: per-record minimizer-prefix quotienting + runtime-selected integer width (uint32/uint64/uint128)
- [X] STRAND-INVARIANT PER-K-MER FRAMING (v0.4.2): exact queries at any m; removed the m>=7 minimizer guard
- [ ] CI/CD WITH GITHUB ACTIONS
- [X] DIFF BETWEEN 2 PRECOMPUTED LISTS (`sskm setop --op diff`; + `diff_size`)
- [X] MERGE OF 2 PRECOMPUTED LISTS (`sskm setop --op union`; + `union_size`)
- [X] INTERSECTION OF 2 PRECOMPUTED LISTS (`sskm setop --op intersection`; + `intersection_size`) — all bucket-parallel (`-t`), output byte-identical; benchmarked vs KMC/CBL/FMSI (scripts/bench/report/SETOPS_REPORT.md)
- [ ] ADD HASHING
- [ ] TEST ON HUMAN, METAGENOME, BACTERIA, RICE PANGENOME?
- [X] TEST AGAINST BQF, SSHASH (benchmark vs bqf/sshash/sbwt/cbl — see scripts/bench)


#### OLD TODOs, TO CHECK
* Add the Bzip2 dependancy in the CMake
* In Skmer header, Check how "pair& operator|= (const kuint& value)" behaves related to "pair(kuint& single) : m_value(0, single)"
* when adding a single kuint, should go to least significant not most.