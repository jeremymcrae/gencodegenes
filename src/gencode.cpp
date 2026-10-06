
#include <algorithm>
#include <cstdint>
#include <deque>
#include <string>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <iostream>

#include "gtf.h"
#include "tx.h"
#include "gencode.h"

namespace gencode {

// check which exon is first, by start position
static bool compareExons(const std::vector<int> & e1, const std::vector<int> & e2) {
    return (e1[0] < e2[0]);
}

static void sort_exons(std::vector<std::vector<int> > & exons) {
    std::sort(exons.begin(), exons.end(), compareExons);
}

// find the index of the exon containing a given chromosome position
static std::uint32_t get_exon_num(const std::vector<std::vector<int> > & exons, int pos) {
    for (std::uint32_t i=0; i<exons.size(); i++) {
        if ((pos >= exons[i][0]) && (pos <= exons[i][1])) {
            return i;
        }
    }
    throw std::invalid_argument("can't find exon for an end codon");
}

// adjust CDS coords for any end codon positions
// 
// the CDS might exclude the start codon and end codon from the CDS range
// we have collected the start and end codon ends in cds_range. We can't
// just set the first CDS coord and last CDS coord to their values though,
// as at least one stop codon spans an intron boundary, which messes up the
// CDS if included as is.
static void include_end_codons(const std::map<std::string, int> & cds_range, TxInfo & info) {
    if (info.cds.size() == 0) {
        return;
    }
    sort_exons(info.cds);
    sort_exons(info.exons);
    int cds_min = cds_range.at("min");
    int cds_max = cds_range.at("max");

    // handle left (5') boundary
    std::uint32_t first_idx = get_exon_num(info.exons, info.cds[0][0]);
    std::uint32_t min_idx = get_exon_num(info.exons, cds_min);
    if (min_idx == first_idx) {
        info.cds[0][0] = cds_min;
    } else {
        info.cds[0][0] = info.exons[first_idx][0];  // extend existing CDS
        std::vector<int> extra_cds = {cds_min, info.exons[min_idx][1]};
        auto it = info.cds.begin();
        info.cds.insert(it, extra_cds);
    }

    // handle right (3') boundary
    std::uint32_t last_idx = get_exon_num(info.exons, info.cds.back()[1]);
    std::uint32_t max_idx = get_exon_num(info.exons, cds_max);
    if (max_idx == last_idx) {
        info.cds.back()[1] = cds_max;
    } else {
        info.cds.back()[1] = info.exons[last_idx][1];  // extend existing CDS
        std::vector<int> extra_cds = {info.exons[max_idx][0], cds_max};
        info.cds.push_back(extra_cds);
    }
}

// set the transcript span from its exons and CDS, for GTFs without transcript lines
static void set_span(TxInfo & info) {
    if (info.start != 0 || info.end != 0) {
        return;
    }
    bool first = true;
    for (auto regions : {&info.exons, &info.cds}) {
        for (auto & x : *regions) {
            info.start = first ? x[0] : std::min(info.start, x[0]);
            info.end = first ? x[1] : std::max(info.end, x[1]);
            first = false;
        }
    }
}

// construct a Tx from the features collected for a transcript
//
// Transcripts with inconsistent coordinates (e.g. a stop codon outside the
// exons) are skipped with a warning, so one malformed transcript doesn't stop
// the rest of the GTF from loading.
static void add_transcript(std::vector<NamedTx> & transcripts, TxInfo & info,
        std::map<std::string, int> & cds_range, std::string & symbol,
        std::vector<std::string> & alt_ids) {
    try {
        // adjust CDS for start and stop codon coords
        include_end_codons(cds_range, info);
        set_span(info);
        Tx tx = Tx(info.name, info.chrom, info.start, info.end, info.strand[0],
            info.transcript_type, info.attributes);
        tx.set_exons(info.exons);
        tx.set_cds(info.cds);
        transcripts.push_back({std::move(symbol), std::move(alt_ids), std::move(tx), info.is_canonical});
    } catch (const std::invalid_argument & e) {
        std::cerr << "skipping transcript " << info.name << ": " << e.what() << std::endl;
    }
}

// features collected for a transcript, while loading the GTF
struct TxEntry {
    TxInfo info;
    std::map<std::string, int> cds_range = {{"max", 0}, {"min", 999999999}};
    std::string symbol;
    std::vector<std::string> alt_ids;
};

// build transcripts in the order they first appeared, freeing each entry once
// its transcript is built
//
// @param pos start of the current GTF line. Transcripts are only built once this
//     passes their end, as no further lines for them can follow in a
//     position-sorted GTF. Use -1 to build all transcripts.
static void build_transcripts(std::vector<NamedTx> & transcripts,
        std::deque<TxEntry> & entries, std::unordered_map<std::string, TxEntry *> & index,
        int pos=-1) {
    while (!entries.empty()) {
        TxEntry & x = entries.front();
        if (pos != -1 && (x.info.end == 0 || pos <= x.info.end)) {
            break;
        }
        index.erase(x.info.name);
        add_transcript(transcripts, x.info, x.cds_range, x.symbol, x.alt_ids);
        entries.pop_front();
    }
}

// collect all features for a transcript into a single object
//
// When we load lines from gencode GTF files, each line represents a single exon
// or CDS, and we need to combine these based on transcript ID. Lines for a
// transcript are usually contiguous, but position-sorted GTFs interleave
// transcripts, so features are collected by transcript ID until the chromosome
// changes (GTFs are grouped by chromosome), or the GTF moves past the transcript.
static void load_transcripts(std::vector<NamedTx> & transcripts, GTF &gtf_file, bool coding=true) {
    std::set<std::string> permit = {"exon", "CDS", "UTR", "transcript", 
        "stop_codon", "start_codon"};
    std::deque<TxEntry> entries;
    std::unordered_map<std::string, TxEntry *> index;
    TxEntry * entry = nullptr;

    GTFLine gtf;

    while (gtf_file.next(gtf)) {
        if (permit.count(gtf.feature) == 0) {
            continue;
        } else if (coding && (gtf.transcript_type != "protein_coding")) {
            continue;
        }

        // only look up the transcript when it differs from the previous line's
        if (entry == nullptr || gtf.tx_id != entry->info.name || gtf.chrom != entry->info.chrom) {
            int pos = (entry != nullptr && gtf.chrom != entry->info.chrom) ? -1 : gtf.start;
            build_transcripts(transcripts, entries, index, pos);
            auto it = index.find(gtf.tx_id);
            if (it != index.end()) {
                entry = it->second;
            } else {
                entries.emplace_back();
                entry = &entries.back();
                index[gtf.tx_id] = entry;
                
                TxInfo & info = entry->info;
                info.name = gtf.tx_id;
                info.chrom = gtf.chrom;
                info.strand = gtf.strand;
                info.is_canonical = gtf.is_canonical;
                info.transcript_type = gtf.transcript_type;
                entry->symbol = gtf.symbol;
                entry->alt_ids = gtf.alternate_ids;
                if (gtf.feature != "transcript") {
                    // without a transcript line, use the first line's attributes,
                    // minus the fields specific to that feature
                    info.attributes = gtf.attributes;
                    for (auto field : {"exon_number", "exon_id", "exon_version"}) {
                        info.attributes.erase(field);
                    }
                }
            }
        }

        TxInfo & info = entry->info;
        std::map<std::string, int> & cds_range = entry->cds_range;
        if (gtf.feature == "transcript") {
            info.start = gtf.start;
            info.end = gtf.end;
            info.attributes = std::move(gtf.attributes);
        } else if (gtf.feature == "CDS") {
            info.cds.push_back(std::vector<int> {gtf.start, gtf.end});
            cds_range["max"] = std::max(std::max(cds_range["max"], gtf.start), gtf.end);
            cds_range["min"] = std::min(std::min(cds_range["min"], gtf.start), gtf.end);
        } else if (gtf.feature == "exon") {
            info.exons.push_back(std::vector<int> {gtf.start, gtf.end});
        } else if ((gtf.feature == "stop_codon") || (gtf.feature == "start_codon")) {
            cds_range["max"] = std::max(std::max(cds_range["max"], gtf.start), gtf.end);
            cds_range["min"] = std::min(std::min(cds_range["min"], gtf.start), gtf.end);
        }
    }

    build_transcripts(transcripts, entries, index);
}

std::vector<NamedTx> open_gencode(std::string path, bool coding) {
    GTF gtf_file(path);
    std::vector<NamedTx> transcripts;
    transcripts.reserve(80000);
    load_transcripts(transcripts, gtf_file, coding);
    return transcripts;
}

bool CompFunc(const GenePoint &l, const GenePoint &r) {
    return l.pos < r.pos;
}

std::vector<std::string> _in_region(std::string chrom, int start, int end, 
        std::map<std::string, std::vector<GenePoint>> & starts, 
        std::map<std::string, std::vector<GenePoint>> & ends,
        int max_window=2500000) {
    
    if (starts.count(chrom) == 0) {
        throw std::invalid_argument("unknown_chrom: " + chrom);
    }
    
    std::vector<GenePoint> & chrom_starts = starts[chrom];
    std::vector<GenePoint> & chrom_ends = ends[chrom];
    std::set<std::size_t> inside;
    std::vector<std::string> symbols;
    symbols.reserve(std::max((end - start) / 50000, 1));  // expect 1 gene / 50 kb 
    
    // find indices to genes with a start inside the region
    int left_idx, right_idx;
    int idx = std::lower_bound(chrom_starts.begin(), chrom_starts.end(), GenePoint {start, "A"}, CompFunc) - chrom_starts.begin();
    left_idx = idx - 1;
    while (idx < (int) chrom_starts.size()) {
        GenePoint & edge = chrom_starts[idx];
        if (edge.pos > end) {
            break;
        } else {
            inside.insert(std::hash<std::string>{}(edge.symbol));
            symbols.push_back(edge.symbol);
        }
        idx += 1;
    }

    // find indices to genes with a end inside the region
    idx = (std::upper_bound(chrom_ends.begin(), chrom_ends.end(), GenePoint {end, "A"}, CompFunc) - chrom_ends.begin());
    right_idx = idx;
    idx = std::min(idx - 1, (int) chrom_ends.size() - 1);
    while (idx >= 0) {
        GenePoint & edge = chrom_ends[idx];
        if (edge.pos < start) {
            break;
        } else if (inside.count(std::hash<std::string>{}(edge.symbol)) == 0) {
            symbols.push_back(edge.symbol);
        }
        idx -= 1;
    }

    if (abs(end - start) > max_window) {
        // if the window is too wide to permit a gene to span it, just return
        return symbols;
    }
    // for genes that encapsulate the region, first find genes that start upstream
    static std::set<std::size_t> starts_before;
    starts_before.clear();
    for (; left_idx>=0; left_idx--) {
        GenePoint & edge = chrom_starts[left_idx];
        if (abs(edge.pos - end) > max_window) { // halt if distant from the region
            break;
        }
        starts_before.insert(std::hash<std::string>{}(edge.symbol));
    }

    // find genes that end downstream of the gene
    int length = (int) chrom_ends.size();
    for (; right_idx<length; right_idx++) {
        GenePoint & edge = chrom_ends[right_idx];
        if (abs(edge.pos - start) > max_window) { // halt if distant from the region
            break;
        }
        if (starts_before.count(std::hash<std::string>{}(edge.symbol)) != 0) {
            symbols.push_back(edge.symbol);
        }
    }

    return symbols;
}

} // namespace

// int main() {
//     std::string path = "/illumina/scratch/deep_learning/public_data/refdata/hg38/genes/gencode.v24.annotation.gtf";
//     gencode::open_gencode(path);
// }
// 
// g++ -std=c++11 gencode.cpp gtf.cpp tx.cpp -lz



