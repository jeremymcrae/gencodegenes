#ifndef GENCODEGENES_GENCODE_H_
#define GENCODEGENES_GENCODE_H_

#include <string>
#include <utility>
#include <vector>

#include "tx.h"

namespace gencode {

// stores HGNC symbol with the transcript, so we can collect transcripts by gene
struct NamedTx {
    std::string symbol;
    std::vector<std::string> alternate_ids;
    Tx tx;
    int canonical_score;
};

// gene coordinates, for finding genes by position
struct GeneInterval {
    int start;
    int end;
    int max_end;      // highest end of this and all preceding intervals
    int idx;          // index of the gene
    int max_end_idx;  // index of the gene ending at max_end
};

// open the gencode annotations GTF, and return Tx objects for each transcript
std::vector<NamedTx> open_gencode(std::string path, bool coding=true);

void index_intervals(std::vector<GeneInterval> & intervals);
std::vector<int> genes_in_region(const std::vector<GeneInterval> & intervals,
    int start, int end);
std::pair<int, int> flanking_genes(const std::vector<GeneInterval> & intervals,
    int pos);

} // namespace

#endif // GENCODEGENES_GENCODE_H_
