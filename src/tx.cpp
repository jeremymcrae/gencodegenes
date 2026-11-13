#include <array>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <stdexcept>
#include <utility>

#include "tx.h"
#include "bases.h"

namespace gencode {

// standard genetic code, indexed by codon with bases encoded as A=0, C=1, G=2,
// T=3, i.e. AAA, AAC, AAG, AAT, ACA, ... TTT
static const char aa_code[] =
    "KNKNTTTTRSRSIIMI"
    "QHQHPPPPRRRRLLLL"
    "EDEDAAAAGGGGVVVV"
    "*Y*YSSSS*CWCLFLF";

// lookup table of 2-bit codes for (upper case) bases, with -1 for other characters
static std::array<int8_t, 256> make_base_codes() {
    std::array<int8_t, 256> table;
    table.fill(-1);
    table['A'] = 0;
    table['C'] = 1;
    table['G'] = 2;
    table['T'] = 3;
    return table;
}

// translate a codon, or return 0 if the codon contains other characters
static char translate_bases(const char * codon) {
    static const std::array<int8_t, 256> codes = make_base_codes();
    int a = codes[(uint8_t) codon[0]];
    int b = codes[(uint8_t) codon[1]];
    int c = codes[(uint8_t) codon[2]];
    if ((a | b | c) < 0) {
        return 0;
    }
    return aa_code[(a << 4) | (b << 2) | c];
}

Tx::Tx(std::string transcript_id, std::string chromosome,
    int start_pos, int end_pos, char strand, std::string _transcript_type,
    std::map<std::string, std::string> _attributes)
    : name(std::move(transcript_id)), chrom(std::move(chromosome)),
      tx_start(start_pos), tx_end(end_pos), tx_strand(strand),
      transcript_type(std::move(_transcript_type)),
      attributes(std::make_shared<std::map<std::string, std::string>>(std::move(_attributes))) {
    if (strand != '+' && strand != '-') {
        throw std::invalid_argument("unknown strand type");
    }
}

// sort order for regions, by start then end
static bool compare_regions(const Region & a, const Region & b) {
    return a.start < b.start || (a.start == b.start && a.end < b.end);
}

// convert [start, end] pairs (e.g. from python) to Regions
static std::vector<Region> to_regions(const std::vector<std::vector<int>> & ranges) {
    std::vector<Region> regions;
    regions.reserve(ranges.size());
    for (auto & range : ranges) {
        if (range.size() != 2) {
            throw std::invalid_argument("regions need a start and end, but got "
                + std::to_string(range.size()) + " values");
        }
        regions.push_back(Region {range[0], range[1]});
    }
    return regions;
}

// total number of bases in a set of regions
static int total_length(const std::vector<Region> & regions) {
    int length = 0;
    for (auto & region : regions) {
        length += (region.end - region.start) + 1;
    }
    return length;
}

static void check_regions(const std::vector<Region> & regions, const std::string & name) {
    for (auto & region : regions) {
        if (region.start > region.end) {
            throw std::invalid_argument(name + " has a region starting after it ends: "
                + std::to_string(region.start) + "-" + std::to_string(region.end));
        }
    }
}

// set exon ranges to the class object
//
// @param exon_ranges list of lists e.g. [[5, 10], [20, 30]]
void Tx::set_exons(std::vector<std::vector<int>> exon_ranges) {
    set_exons(to_regions(exon_ranges));
}

void Tx::set_exons(std::vector<Region> exon_ranges) {
    check_regions(exon_ranges, get_name());
    std::sort(exon_ranges.begin(), exon_ranges.end(), compare_regions);
    exons = std::move(exon_ranges);
}

// set CDS ranges to the class object
//
// @param cds_ranges nested vector of ints e.g. [[5, 10], [20, 30]]
void Tx::set_cds(std::vector<std::vector<int>> cds_ranges) {
    set_cds(to_regions(cds_ranges));
}

void Tx::set_cds(std::vector<Region> cds_ranges) {
    check_regions(cds_ranges, get_name());
    cds.clear();
    cds_sequence.clear();
    cds_length = 0;
    cds_min = 0;
    cds_max = 0;
    if (cds_ranges.empty()) {
        return;
    }
    
    // without exons, a single CDS region within the transcript gets an exon
    // spanning the transcript
    if (exons.empty()) {
        Region region = cds_ranges[0];
        if (cds_ranges.size() > 1 || region.start < tx_start || region.end > tx_end) {
            throw std::invalid_argument(get_name() + " lacks exon coordinates");
        }
        exons.push_back(Region {tx_start, tx_end});
    }
    
    std::sort(cds_ranges.begin(), cds_ranges.end(), compare_regions);
    cds = std::move(cds_ranges);
    cds_min = cds.front().start;
    cds_max = cds.back().end;
    move_cds_ends_into_exons();
    cds_length = total_length(cds);
    
    if (!genomic_sequence.empty()) {
        // rebuild the CDS sequence, which can also complete the final codon
        add_genomic_sequence(genomic_sequence, true);
    }
}

// The CDS coordinates are repaired in two steps: CDS ends lying in introns are
// moved into the exons when the CDS is set, then an incomplete final codon is
// completed once the genomic sequence is added (as that needs the bases past the
// CDS end). The start and stop codons are also included when loading the GTF.

// fix CDS start and ends which occur outside exons.
//
// The main cause is when the stop codon overlaps an exon boundary. In 
// some of these cases, the CDS end is placed within the adjacent intron 
// sequence (if that allows stop site sequence). To fix this, we adjust the 
// position to the offset distance within the adjacent exon.
//
// @param position chromosome position
//
// @returns adjusted chromosome positions
Region Tx::fix_cds_boundary(int position) {
    if ( is_exonic(position) ) {
        throw std::invalid_argument( "position shouldn't be in exons" );
    }
    
    Region region = get_closest_exon(position);
    
    int start_dist = std::abs(position - region.start);
    int end_dist = std::abs(position - region.end);
    
    // if we are closer to the start, then we go back an exon, otherwise forward
    int i = (start_dist < end_dist) ? closest_exon_num(region.start) - 1
        : closest_exon_num(region.end) + 1;
    if (i < 0 || i >= (int) exons.size()) {
        throw std::invalid_argument(get_name() + " CDS lies outside the exons");
    }
    
    // place the bases which ran into the intron at the facing end of that exon
    int start;
    int end;
    if (start_dist < end_dist) {
        end = exons[i].end;
        start = end - start_dist;
    } else {
        start = exons[i].start;
        end = start + end_dist;
    }
    
    return Region {start, end};
}

// move CDS ends which lie in introns into the adjacent exons. The intronic part
// is trimmed from the end CDS region, and added as a new region in the adjacent
// exon (see fix_cds_boundary).
void Tx::move_cds_ends_into_exons() {
    if (!is_exonic(cds_min)) {
        Region region = fix_cds_boundary(cds_min);
        cds.front().start += region.end - region.start;
        cds.insert(cds.begin(), region);
        cds_min = region.start;
    }
    if (!is_exonic(cds_max)) {
        Region region = fix_cds_boundary(cds_max);
        cds.back().end -= region.end - region.start;
        cds.push_back(region);
        cds_max = region.end;
    }
}

// extend the CDS by 1-2 bases to complete the final codon, if it is incomplete.
//
// Some transcripts terminate in the middle of a codon, which causes problems
// when translating the final codon. This is rare, and it is rarer that we
// select these transcripts for analysis. TNS3 is an example, where the de novos
// in the gene can only be contained within a transcript that is incomplete at
// the 3' end. The exons are also extended if the CDS now runs past them.
void Tx::complete_final_codon() {
    int missing = (3 - cds_sequence.size() % 3) % 3;
    if (cds.empty() || missing == 0) {
        return;
    }
    
    // the genomic sequence is on the + strand, so the extra bases come from
    // just past the 3' CDS end on the + strand, or just before it (then
    // reverse complemented) on the - strand
    bool fwd = is_forward();
    int end = get_cds_end();
    int idx = genomic_index(fwd ? end + 1 : end - missing);
    if (idx < 0 || idx + missing > (int) genomic_sequence.size()) {
        // the genomic sequence doesn't extend far enough to complete the codon
        return;
    }
    std::string extra = genomic_sequence.substr(idx, missing);
    if (fwd) {
        cds.back().end += missing;
        cds_sequence += extra;
    } else {
        cds.front().start -= missing;
        cds_sequence += reverse_complement(extra);
    }
    cds_min = cds.front().start;
    cds_max = cds.back().end;
    cds_length = total_length(cds);
    
    exons.front().start = std::min(exons.front().start, cds_min);
    exons.back().end = std::max(exons.back().end, cds_max);
}

int Tx::get_cds_start() {
    return is_forward() ? cds_min : cds_max;
}

int Tx::get_cds_end() {
    return is_forward() ? cds_max : cds_min;
}

// checks if a position lies within the exon ranges
//
// @param position integer chromosome position e.g. 10000000
bool Tx::is_exonic(int pos) {
    Region exon = get_closest_exon(pos);
    return is_exonic(pos, exon);
}

// check if a position lies within an exon
//
// @param position integer chromosome position e.g. 10000000
// @param exon Region for an exon
bool Tx::is_exonic(int pos, Region exon) {
    return (pos >= exon.start) && (pos <= exon.end);
}

// check if a region starts before a position, for binary searches of sorted regions
static bool starts_before(const Region & region, int pos) {
    return region.start < pos;
}

// find the index of the closest exon to a chromosome position
//
// @param pos chromosomal nucleotide position
// @returns index of the exon containing the position, or the nearest exon
int Tx::closest_exon_num(int pos) {
    size_t idx = std::lower_bound(exons.begin(), exons.end(), pos, starts_before) - exons.begin();
    if (idx == 0) {
        return idx;
    } else if (idx == exons.size()) {
        return idx - 1;
    }
    
    // get the exon, and the preceding exon (lower bound finds exon after pos)
    Region a = exons[idx - 1];
    Region b = exons[idx];
    
    // if the position is within an exon, return that exon number
    if ((pos >= a.start) && (pos <= a.end)) {
        return idx - 1;
    } else if ((pos >= b.start) && (pos <= b.end)) {
        return idx;
    }
    
    // pos lies in the intron between a and b, so pick the closer exon
    int delta_a = std::min(std::abs(pos - a.start), std::abs(pos - a.end));
    int delta_b = std::min(std::abs(pos - b.start), std::abs(pos - b.end));
    
    if (delta_a < delta_b) {
        return idx - 1;
    }
    
    return idx;
}

// find the closest exon for a chromosome position
//
// @param position integer chromosome position e.g. 10000000
Region Tx::get_closest_exon(int pos) {
    if (exons.size() == 0) {
        throw std::invalid_argument("no exons assigned to this transcript");
    }
    int idx = closest_exon_num(pos);
    return exons[idx];
}

// checks if a position lies within the CDS
//
// @param position integer chromosome position e.g. 10000000
bool Tx::in_coding_region(int pos) {
    return overlaps_cds(pos, pos);
}

// shift an intronic site to the nearest exon boundary (but track the offset)
CDS_coords Tx::to_closest_exon(int pos, Region exon) {
    int site = (std::abs(exon.start - pos) < std::abs(exon.end - pos)) ? exon.start : exon.end;
    int offset = is_forward() ? pos - site : site - pos;
    return CDS_coords {site, offset};
}

// count the exonic bases before a position, on the + strand
//
// @param pos chromosome position
// @param idx index of the exon to measure the position from
int Tx::exonic_position(int pos, int idx) {
    int bases = pos - exons[idx].start;
    for (int i = 0; i < idx; i++) {
        bases += (exons[i].end - exons[i].start) + 1;
    }
    return bases;
}

// find the distance in exonic bp from the CDS start to a position, in the
// direction of transcription. Intronic positions are measured from the
// closest boundary of the closest exon, with the distance into the intron
// given as the offset.
//
// @param pos chromosome position
// @returns CDS_coords with the exonic distance and the intronic offset
CDS_coords Tx::get_coding_distance(int pos) {
    if (exons.size() == 0) {
        throw std::invalid_argument("no exons assigned to this transcript");
    }
    int idx = closest_exon_num(pos);
    int offset = 0;
    if (!is_exonic(pos, exons[idx])) {
        CDS_coords site = to_closest_exon(pos, exons[idx]);
        offset = site.offset;
        pos = site.position;
    }
    
    // count exonic bases from the CDS start, then flip the sign on the - strand,
    // where transcription runs towards lower positions
    int cds_start = get_cds_start();
    int delta = exonic_position(pos, idx)
        - exonic_position(cds_start, closest_exon_num(cds_start));
    return CDS_coords {is_forward() ? delta : -delta, offset};
}

// find the chromosome position of a CDS position
//
// @param cds_position distance in coding bp from the CDS start
// @param offset distance into the adjacent intron, in the direction of
//     transcription
// @returns chromosome position
int Tx::get_position_on_chrom(int cds_position, int offset) {
    // walk the CDS regions in the direction of transcription, until reaching
    // the region containing the position
    bool fwd = is_forward();
    int remaining = cds_position;
    for (size_t i = 0; i < cds.size(); i++) {
        const Region & region = fwd ? cds[i] : cds[cds.size() - 1 - i];
        int length = (region.end - region.start) + 1;
        if (remaining >= 0 && remaining < length) {
            return fwd ? region.start + remaining + offset : region.end - remaining - offset;
        }
        remaining -= length;
    }
    throw std::invalid_argument("position not in CDS");
}

// figure out the codon position for a position
int Tx::get_codon_number_for_cds_position(int cds_position) {
    // floor division, so negative (upstream) positions round down
    return (cds_position >= 0) ? cds_position / 3 : (cds_position - 2) / 3;
}

// get the position within a codon (in 0 based format eg 0, 1, 2)
int Tx::get_position_within_codon(int cds_position) {
    // matches the floor division in get_codon_number_for_cds_position
    return ((cds_position % 3) + 3) % 3;
}

// add the CDS sequence
void Tx::add_cds_sequence(std::string cds_dna) {
    to_upper(cds_dna);
    cds_sequence = std::move(cds_dna);
}

// set how far the genomic sequence extends beyond each end of the transcript.
// This can't change once a sequence is loaded, as it would shift every position.
void Tx::set_genomic_offset(int offset) {
    if (!genomic_sequence.empty() && offset != gdna_offset) {
        throw std::invalid_argument("cannot change the genomic offset of " + get_name()
            + " from " + std::to_string(gdna_offset) + " to " + std::to_string(offset)
            + " after its genomic sequence is loaded");
    }
    gdna_offset = offset;
}

// add the genomic sequence, and build the CDS sequence from it. This completes
// an incomplete final codon, so can extend the CDS and exon coordinates.
//
// @param gdna string for genomic sequence. If the transcript strand is '-',
//    then the DNA sequence will be for the - strand, so we need to
//    reorient the DNA to the + strand.
// @param plus_strand true if gdna is already on the + strand
void Tx::add_genomic_sequence(std::string gdna, bool plus_strand) {
    int expected = (tx_end - tx_start) + 1 + 2 * gdna_offset;
    if ((int) gdna.size() != expected) {
        throw std::invalid_argument(get_name() + " genomic sequence is "
            + std::to_string(gdna.size()) + " bp, but should be " + std::to_string(expected)
            + " bp (the transcript span, plus the offset on each side)");
    }
    to_upper(gdna);
    if (!is_forward() && !plus_strand) {
        // orient the DNA sequence to the + strand.
        gdna = reverse_complement(gdna);
    }
    
    genomic_sequence = std::move(gdna);
    extract_cds_sequence();
    complete_final_codon();
}

// build the CDS sequence from the genomic sequence, checking it matches any CDS
// sequence set already
void Tx::extract_cds_sequence() {
    std::string cds_seq;
    for (auto &region : cds) {
        int x = genomic_index(region.start);
        int len = (region.end - region.start) + 1;
        cds_seq += genomic_sequence.substr(x, len);
    }
    if (!is_forward()) {
        cds_seq = reverse_complement(cds_seq);
    }
    
    // don't check the CDS matches expectations if the CDS sequence isn't set
    if (cds_sequence == "") {
        cds_sequence = cds_seq;
    } else if (cds_seq != cds_sequence) {
        // do a sanity check to check that we've got the right cds sequence, this
        // fails for at least one gene (CCDC18), which begins with a N, and
        // throws the coordinates off
        std::string msg = "Coding sequence from gene coordinates doesn't match "
            "coding sequence obtained from Ensembl.\nTranscript:" + get_name() +
            "\n" + cds_seq + "\n\nshould be\n" + cds_sequence + "\n";
        throw std::invalid_argument(msg);
    }
}

// reverse complement a DNA or RNA sequence
std::string Tx::reverse_complement(const std::string & seq) {
    size_t n = seq.size();
    std::string revcomp(n, 'N');
    for (size_t i = 0; i < n; i++) {
        revcomp[i] = complement(seq[n - 1 - i]);
    }
    return revcomp;
}

// obtains the sequence around a position
//
// @param pos integer for chromosomal nucleotide position.
// @param length integer for length of region to select.
//
// @returns nucleotide sequence e.g. a trinucleotide such as 'TCC',
//    centered around the required position, given with respect to the fwd
//    strand.
std::string Tx::get_centered_sequence(int pos, int length) {
    if (pos < 0) {
        throw std::invalid_argument( "sequence position < 0" );
    }
    if (length % 2 != 1) {
        throw std::invalid_argument( "length is not an odd number") ;
    }
    
    int start = pos - length / 2;
    return get_seq_in_region(start, start + length - 1);
}

// get the + strand genomic sequence from start to end
std::string Tx::get_seq_in_region(int start, int end) {
    if (genomic_sequence.empty()) {
        throw std::invalid_argument(get_name() + " lacks genomic sequence (e.g. "
            "from a genome fasta)");
    }
    if (start < (get_start() - gdna_offset) || end > (get_end() + gdna_offset)) {
        throw std::invalid_argument( "sequence position not in gene range" );
    }
    return genomic_sequence.substr(genomic_index(start), (end - start) + 1);
}

// get the codon sequence for a given codon_number
std::string Tx::get_codon_sequence(int codon) { 
    if (cds_sequence.empty()) {
        throw std::invalid_argument(get_name() + " lacks coding sequence (e.g. "
            "from a genome fasta)");
    }
    if (codon < 0) {
        throw std::invalid_argument( "codon position < 0" );
    }
    if (codon >= (int)cds_sequence.size() / 3) {
        throw std::invalid_argument( "codon position not in gene range" );
    }
    
    return cds_sequence.substr(codon * 3, 3);
}

// translate a single codon to amino acid
// @param seq codon sequence
// @returns translated amino acid sequence as single char
char Tx::translate_codon(const std::string & codon) {
    if (codon.size() != 3) {
        throw std::invalid_argument("translate_codon() needs a 3 bp sequence");
    }
    char aa = translate_bases(codon.data());
    if (aa == 0) {
        throw std::invalid_argument("cannot translate codon: " + codon);
    }
    return aa;
}

// translate a DNA codon to a single character amino acid
//
// @param seq codon sequence, or longer DNA sequences.
// @returns translated amino acid sequence as single character codes.
std::string Tx::translate(std::string seq) {
    to_upper(seq);
    if (seq.size() % 3 != 0) {
        throw std::invalid_argument("translate() needs a multiple of 3 bp, but got "
            + std::to_string(seq.size()) + " bp");
    }
    
    std::string protein(seq.size() / 3, '\0');
    for (size_t i = 0; i < protein.size(); i++) {
        protein[i] = translate_bases(seq.data() + i * 3);
        if (protein[i] == 0) {
            throw std::invalid_argument("cannot translate codon: " + seq.substr(i * 3, 3));
        }
    }
    return protein;
}

// get the details of the codon which a variant resides in
//
// @param bp nucleotide position of the variant (within the transcript gene range)
// @returns dictionary of codon sequence, cds position, amino acid that the
//    codon translates to, and position within the codon.
Codon Tx::get_codon_info(int bp) {
    bool in_coding = in_coding_region(bp);
    CDS_coords site = get_coding_distance(bp);
    
    if ((site.position < 0) || (site.position >= cds_length)) {
        throw std::invalid_argument( "position not inside CDS region" );
    }
    
    // intronic sites between the CDS start and end have a CDS position, but
    // aren't in a codon
    Codon codon {site.position, "", -1, -1, '\0', site.offset, in_coding};
    if (in_coding) {
        codon.codon_number = get_codon_number_for_cds_position(site.position);
        codon.intra_codon = get_position_within_codon(site.position);
        codon.codon_seq = get_codon_sequence(codon.codon_number);
        codon.initial_aa = translate_codon(codon.codon_seq);
    }
    return codon;
}

// get the distance in bp for a variant to the nearest exon boundary
//
// @param bp nucleotide position of the variant (within the transcript gene range)
// @returns distance in base-pairs to the nearest exon boundary.
int Tx::get_boundary_distance(int bp) {
    Region exon = get_closest_exon(bp);
    int distance = std::min(std::abs(exon.start - bp), std::abs(exon.end - bp));
    
    // sites within the coding region are actually one bp further away,
    // since we are measuring the distance to the base inside the exon
    // boundary
    if ( in_coding_region(bp) ) {
        distance += 1;
    }
    
    return distance;
}

// find the index of the first region overlapping a genomic range. Returns -1
// if none overlap.
int Tx::overlapping_region_num(int start, int end, std::vector<Region> & group) {
    // find the first region ending at or after the start
    auto it = std::lower_bound(group.begin(), group.end(), start,
        [](const Region & r, int pos) { return r.end < pos; });
    if (it == group.end() || !overlap(start, end, it->start, it->end)) {
        return -1;
    }
    return (int) (it - group.begin());
}

// checks if a genomic range overlaps any CDS range
bool Tx::overlaps_cds(int start, int end) {
    return overlapping_region_num(start, end, cds) >= 0;
}

} // namespace
