
#include <algorithm>
#include <climits>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "gtf.h"
#include "tx.h"
#include "gencode.h"

namespace gencode {

// features collected for a transcript while loading the GTF, to construct a Tx
struct TxInfo {
    std::string name;
    std::string chrom;
    int start = 0;
    int end = 0;
    char strand;
    std::string transcript_type;
    std::vector<Region> exons;
    std::vector<Region> cds;
    // CDS range, including the start and stop codons
    int cds_min = INT_MAX;
    int cds_max = 0;
    int is_canonical = 0;
    std::string symbol;
    std::vector<std::string> alt_ids;
    std::map<std::string, std::string> attributes;
};

// find the index of the exon containing a given chromosome position
static size_t get_exon_num(const std::vector<Region> & exons, int pos) {
    for (size_t i = 0; i < exons.size(); i++) {
        if ((pos >= exons[i].start) && (pos <= exons[i].end)) {
            return i;
        }
    }
    throw std::invalid_argument("can't find exon for an end codon");
}

// adjust CDS coords for any end codon positions
// 
// the CDS might exclude the start codon and end codon from the CDS range
// we have collected the start and end codon ends in cds_min and cds_max. We
// can't just set the first CDS coord and last CDS coord to their values though,
// as at least one stop codon spans an intron boundary, which messes up the
// CDS if included as is.
static void include_end_codons(TxInfo & info) {
    if (info.cds.empty()) {
        return;
    }
    auto by_start = [](const Region & a, const Region & b) { return a.start < b.start; };
    std::sort(info.cds.begin(), info.cds.end(), by_start);
    std::sort(info.exons.begin(), info.exons.end(), by_start);

    // handle left (5') boundary
    size_t first_idx = get_exon_num(info.exons, info.cds[0].start);
    size_t min_idx = get_exon_num(info.exons, info.cds_min);
    if (min_idx == first_idx) {
        info.cds[0].start = info.cds_min;
    } else {
        info.cds[0].start = info.exons[first_idx].start;  // extend existing CDS
        info.cds.insert(info.cds.begin(), Region {info.cds_min, info.exons[min_idx].end});
    }

    // handle right (3') boundary
    size_t last_idx = get_exon_num(info.exons, info.cds.back().end);
    size_t max_idx = get_exon_num(info.exons, info.cds_max);
    if (max_idx == last_idx) {
        info.cds.back().end = info.cds_max;
    } else {
        info.cds.back().end = info.exons[last_idx].end;  // extend existing CDS
        info.cds.push_back(Region {info.exons[max_idx].start, info.cds_max});
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
            info.start = first ? x.start : std::min(info.start, x.start);
            info.end = first ? x.end : std::max(info.end, x.end);
            first = false;
        }
    }
}

// construct a Tx from the features collected for a transcript
//
// Transcripts with inconsistent coordinates (e.g. a stop codon outside the
// exons) are skipped with a warning, so one malformed transcript doesn't stop
// the rest of the GTF from loading.
static void add_transcript(std::vector<NamedTx> & transcripts, TxInfo & info) {
    try {
        include_end_codons(info);
        set_span(info);
        Tx tx(info.name, info.chrom, info.start, info.end, info.strand,
            info.transcript_type, std::move(info.attributes));
        tx.set_exons(std::move(info.exons));
        tx.set_cds(std::move(info.cds));
        transcripts.push_back({std::move(info.symbol), std::move(info.alt_ids),
            std::move(tx), info.is_canonical});
    } catch (const std::invalid_argument & e) {
        std::cerr << "skipping transcript " << info.name << ": " << e.what() << std::endl;
    }
}

// build transcripts in the order they first appeared, freeing each entry once
// its transcript is built
//
// @param pos start of the current GTF line. Transcripts are only built once this
//     passes their end, as no further lines for them can follow in a
//     position-sorted GTF. Use -1 to build all transcripts.
static void build_transcripts(std::vector<NamedTx> & transcripts,
        std::deque<TxInfo> & entries, std::unordered_map<std::string, TxInfo *> & index,
        int pos=-1) {
    while (!entries.empty()) {
        TxInfo & info = entries.front();
        if (pos != -1 && (info.end == 0 || pos <= info.end)) {
            break;
        }
        index.erase(info.name);
        add_transcript(transcripts, info);
        entries.pop_front();
    }
}

// find the entry for a GTF line's transcript, or start a new one from the line
static TxInfo * get_entry(std::deque<TxInfo> & entries,
        std::unordered_map<std::string, TxInfo *> & index, GTFLine & gtf) {
    auto it = index.find(gtf.tx_id);
    if (it != index.end()) {
        return it->second;
    }
    entries.emplace_back();
    TxInfo * info = &entries.back();
    index[gtf.tx_id] = info;

    info->name = gtf.tx_id;
    info->chrom = gtf.chrom;
    info->strand = gtf.strand[0];
    info->is_canonical = gtf.is_canonical;
    info->transcript_type = gtf.transcript_type;
    info->symbol = gtf.symbol;
    info->alt_ids = gtf.alternate_ids;
    if (gtf.feature != "transcript") {
        // without a transcript line, use the first line's attributes, minus the
        // fields specific to that feature
        info->attributes = std::move(gtf.attributes);
        for (auto field : {"exon_number", "exon_id", "exon_version"}) {
            info->attributes.erase(field);
        }
    }
    return info;
}

enum class Feature { other, transcript, exon, cds, codon, utr };

static Feature get_feature(const std::string & feature) {
    if (feature == "exon") { return Feature::exon; }
    if (feature == "CDS") { return Feature::cds; }
    if (feature == "UTR") { return Feature::utr; }
    if (feature == "transcript") { return Feature::transcript; }
    if (feature == "start_codon" || feature == "stop_codon") { return Feature::codon; }
    return Feature::other;
}

// collect all features for a transcript into a single object
//
// When we load lines from gencode GTF files, each line represents a single exon
// or CDS, and we need to combine these based on transcript ID. Lines for a
// transcript are usually contiguous, but position-sorted GTFs interleave
// transcripts, so features are collected by transcript ID until the chromosome
// changes (GTFs are grouped by chromosome), or the GTF moves past the transcript.
static void load_transcripts(std::vector<NamedTx> & transcripts, GTF & gtf_file, bool coding) {
    std::deque<TxInfo> entries;
    std::unordered_map<std::string, TxInfo *> index;
    TxInfo * info = nullptr;
    GTFLine gtf;

    while (gtf_file.next(gtf)) {
        Feature feature = get_feature(gtf.feature);
        if (feature == Feature::other || (coding && gtf.transcript_type != "protein_coding")) {
            continue;
        }

        // only look up the transcript when it differs from the previous line's
        if (info == nullptr || gtf.tx_id != info->name || gtf.chrom != info->chrom) {
            bool new_chrom = info != nullptr && gtf.chrom != info->chrom;
            build_transcripts(transcripts, entries, index, new_chrom ? -1 : gtf.start);
            info = get_entry(entries, index, gtf);
        }

        if (feature == Feature::transcript) {
            info->start = gtf.start;
            info->end = gtf.end;
            info->attributes = std::move(gtf.attributes);
        } else if (feature == Feature::exon) {
            info->exons.push_back(Region {gtf.start, gtf.end});
        } else if (feature == Feature::cds) {
            info->cds.push_back(Region {gtf.start, gtf.end});
        }
        if (feature == Feature::cds || feature == Feature::codon) {
            info->cds_min = std::min({info->cds_min, gtf.start, gtf.end});
            info->cds_max = std::max({info->cds_max, gtf.start, gtf.end});
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

// find genes overlapping a region
//
// @param starts gene start positions per chromosome, sorted by position
// @param ends gene end positions per chromosome, sorted by position
// @param max_window maximum gene length, which limits how far beyond the region
//     we search for genes that span the whole region
// @returns keys for the genes overlapping the region
std::vector<std::string> _in_region(std::string chrom, int start, int end, 
        std::map<std::string, std::vector<GenePoint>> & starts, 
        std::map<std::string, std::vector<GenePoint>> & ends,
        int max_window) {
    auto starts_it = starts.find(chrom);
    if (starts_it == starts.end()) {
        throw std::invalid_argument("unknown_chrom: " + chrom);
    }
    const std::vector<GenePoint> & chrom_starts = starts_it->second;
    const std::vector<GenePoint> & chrom_ends = ends.at(chrom);
    const int n_starts = (int) chrom_starts.size();
    const int n_ends = (int) chrom_ends.size();
    std::hash<std::string> hash;
    std::vector<std::string> symbols;

    // find genes with a start inside the region
    std::unordered_set<std::size_t> inside;
    int first_start = std::lower_bound(chrom_starts.begin(), chrom_starts.end(),
        GenePoint {start, ""}, CompFunc) - chrom_starts.begin();
    for (int i = first_start; i < n_starts && chrom_starts[i].pos <= end; i++) {
        inside.insert(hash(chrom_starts[i].symbol));
        symbols.push_back(chrom_starts[i].symbol);
    }

    // find genes with an end inside the region (that don't also start inside it)
    int after_end = std::upper_bound(chrom_ends.begin(), chrom_ends.end(),
        GenePoint {end, ""}, CompFunc) - chrom_ends.begin();
    for (int i = after_end - 1; i >= 0 && chrom_ends[i].pos >= start; i--) {
        if (inside.count(hash(chrom_ends[i].symbol)) == 0) {
            symbols.push_back(chrom_ends[i].symbol);
        }
    }

    if (std::abs(end - start) > max_window) {
        // if the window is too wide to permit a gene to span it, just return
        return symbols;
    }

    // find genes that span the region, i.e. start upstream and end downstream
    std::unordered_set<std::size_t> starts_before;
    for (int i = first_start - 1; i >= 0 && std::abs(chrom_starts[i].pos - end) <= max_window; i--) {
        starts_before.insert(hash(chrom_starts[i].symbol));
    }
    for (int i = after_end; i < n_ends && std::abs(chrom_ends[i].pos - start) <= max_window; i++) {
        if (starts_before.count(hash(chrom_ends[i].symbol)) != 0) {
            symbols.push_back(chrom_ends[i].symbol);
        }
    }

    return symbols;
}

} // namespace
