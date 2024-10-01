#include "../include/dbg/datastruct/bit.hpp"

namespace dbglib::bit {

uint64_t pdep_u64(uint64_t a, uint64_t mask) {
    uint64_t dst = 0;
    std::size_t m = 0;
    while (m < 64) {
        if ((mask & (1 << m)) != 0) {
            dst |= (a & static_cast<uint64_t>(1)) << m;
            a >>= 1;
        }
        ++m;
    }
    return dst;
}

}