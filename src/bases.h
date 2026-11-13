#ifndef GENCODEGENES_BASES_H_
#define GENCODEGENES_BASES_H_

#include <array>
#include <cstdint>
#include <string>

namespace gencode {

// check for overlap of two ranges
inline bool overlap(int start1, int end1, int start2, int end2) {
    return (end1 >= start2) && (end2 >= start1);
}

inline void to_upper(std::string & seq) {
    for (auto & base : seq) {
        base = (base >= 'a' && base <= 'z') ? (char) (base - ('a' - 'A')) : base;
    }
}

// lookup table of complementary bases, including IUPAC ambiguity codes. Any
// other character complements to 'N'.
inline std::array<char, 256> make_complements() {
    std::array<char, 256> table;
    table.fill('N');
    const std::string bases = "ACGTUNRYKMSWBDHVacgtunrykmswbdhv";
    const std::string complements = "TGCAANYRMKSWVHDBtgcaanyrmkswvhdb";
    for (size_t i = 0; i < bases.size(); i++) {
        table[(uint8_t) bases[i]] = complements[i];
    }
    return table;
}

// complement a single base
inline char complement(char base) {
    static const std::array<char, 256> complements = make_complements();
    return complements[(uint8_t) base];
}

} // namespace

#endif  // GENCODEGENES_BASES_H_
