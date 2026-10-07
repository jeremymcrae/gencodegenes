#ifndef GENCODEGENES_GENCODE_H_
#define GENCODEGENES_GENCODE_H_

#include <string>
#include <map>
#include <vector>

#include "tx.h"

namespace gencode {

// stores HGNC symbol with the transcript, so we can collect transcripts by gene
struct NamedTx {
    std::string symbol;
    std::vector<std::string> alternate_ids;
    Tx tx;
    int is_canonical;
};

struct GenePoint {
    int pos;
    std::string symbol;
};

// open the gencode annotations GTF, and return Tx objects for each transcript
std::vector<NamedTx> open_gencode(std::string path, bool coding=true);

bool CompFunc(const GenePoint &l, const GenePoint &r);
std::vector<std::string> _in_region(std::string chrom, int start, int end, 
    std::map<std::string, std::vector<GenePoint>> & starts, 
    std::map<std::string, std::vector<GenePoint>> & ends,
    int max_window);

} // namespace

#endif // GENCODEGENES_GENCODE_H_
