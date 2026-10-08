#include <array>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "tx.h"

namespace gencode {

static std::unordered_map<std::string, char> aa_code = {
    {"AAA", 'K'}, {"AAC", 'N'}, {"AAG", 'K'}, {"AAT", 'N'},
    {"ACA", 'T'}, {"ACC", 'T'}, {"ACG", 'T'}, {"ACT", 'T'},
    {"AGA", 'R'}, {"AGC", 'S'}, {"AGG", 'R'}, {"AGT", 'S'},
    {"ATA", 'I'}, {"ATC", 'I'}, {"ATG", 'M'}, {"ATT", 'I'},
    {"CAA", 'Q'}, {"CAC", 'H'}, {"CAG", 'Q'}, {"CAT", 'H'},
    {"CCA", 'P'}, {"CCC", 'P'}, {"CCG", 'P'}, {"CCT", 'P'},
    {"CGA", 'R'}, {"CGC", 'R'}, {"CGG", 'R'}, {"CGT", 'R'},
    {"CTA", 'L'}, {"CTC", 'L'}, {"CTG", 'L'}, {"CTT", 'L'},
    {"GAA", 'E'}, {"GAC", 'D'}, {"GAG", 'E'}, {"GAT", 'D'},
    {"GCA", 'A'}, {"GCC", 'A'}, {"GCG", 'A'}, {"GCT", 'A'},
    {"GGA", 'G'}, {"GGC", 'G'}, {"GGG", 'G'}, {"GGT", 'G'},
    {"GTA", 'V'}, {"GTC", 'V'}, {"GTG", 'V'}, {"GTT", 'V'},
    {"TAA", '*'}, {"TAC", 'Y'}, {"TAG", '*'}, {"TAT", 'Y'},
    {"TCA", 'S'}, {"TCC", 'S'}, {"TCG", 'S'}, {"TCT", 'S'},
    {"TGA", '*'}, {"TGC", 'C'}, {"TGG", 'W'}, {"TGT", 'C'},
    {"TTA", 'L'}, {"TTC", 'F'}, {"TTG", 'L'}, {"TTT", 'F'}};

// Constructor for Tx class
Tx::Tx(std::string transcript_id, std::string chromosome,
    int start_pos, int end_pos, char strand, std::string _transcript_type,
    std::map<std::string, std::string> _attributes) {
    
    name = std::move(transcript_id);
    
    chrom = std::move(chromosome);
    tx_start = start_pos;
    tx_end = end_pos;
    
    if ( strand != '+' && strand != '-' ) {
        throw std::invalid_argument( "unknown strand type" );
    }
    tx_strand = strand;
    transcript_type = std::move(_transcript_type);
    attributes = std::move(_attributes);
}

static bool compare_regions(const Region & a, const Region & b) {
    return a.start < b.start || (a.start == b.start && a.end < b.end);
}

static std::vector<Region> to_regions(const std::vector<std::vector<int>> & ranges) {
    std::vector<Region> regions;
    regions.reserve(ranges.size());
    for (auto & range : ranges) {
        regions.push_back(Region {range[0], range[1]});
    }
    return regions;
}

// set exon ranges to the class object
//
// @param exon_ranges list of lists e.g. [[5, 10], [20, 30]]
void Tx::set_exons(std::vector<std::vector<int>> exon_ranges) {
    set_exons(to_regions(exon_ranges));
}

void Tx::set_exons(std::vector<Region> exon_ranges) {
    std::sort(exon_ranges.begin(), exon_ranges.end(), compare_regions);
    exons = std::move(exon_ranges);
    exon_to_cds.clear();
}

// set CDS ranges to the class object
//
// @param cds_ranges nested vector of ints e.g. [[5, 10], [20, 30]]
void Tx::set_cds(std::vector<std::vector<int>> cds_ranges) {
    set_cds(to_regions(cds_ranges));
}

void Tx::set_cds(std::vector<Region> cds_ranges) {
    cds.clear();
    exon_to_cds.clear();
    if (cds_ranges.empty()) {
        // no subsequent adjustments for transcripts without CDS
        return;
    }
    
    // If the transcript lacks exon coordinates and only has a single CDS
    // region, then if the CDS region fits within the gene range, make a
    // single exon, using the transcript start and end. This prevents issues
    // when adding the CDS coordinates for the transcript.
    if ( exons.empty() ) {
        Region exon = cds_ranges[0];
        if (cds_ranges.size() == 1 && std::min(exon.start, exon.end) >= tx_start \
                && std::max(exon.start, exon.end) <= tx_end) {
            Region region {tx_start, tx_end};
            exons.push_back(region);
        } else {
            std::string msg = get_name() + " lacks exon coordinates";
            throw std::invalid_argument( msg );
        }
    }
    
    std::sort(cds_ranges.begin(), cds_ranges.end(), compare_regions);
    cds = std::move(cds_ranges);
    
    cds_min = cds[0].start;
    cds_max = cds.back().end;
    
    // sometimes the cds start and end lie outside the exons. We find the
    // offset to the closest exon boundary, and use that to place the
    // position within the up or downstream exon, and correct the initial
    // cds boundaries
    if (!is_exonic(cds_min)) {
        Region region = fix_cds_boundary(cds_min);
        int delta = std::abs(region.end - region.start);
        
        cds[0] = Region {cds[0].start + delta, cds[0].end};
        cds_min = region.start;
        cds.insert(cds.begin(), region);
    }
    
    if (!is_exonic(cds_max)) {
        Region region = fix_cds_boundary(cds_max);
        int delta = std::abs(region.end - region.start);
        
        size_t idx = cds.size() - 1;
        cds[idx] = Region {cds[idx].start, cds[idx].end - delta};
        cds_max = region.end;
        cds.push_back(region);
    }
    
    cds_length = 0;
    for (auto region : cds) {
        cds_length += (region.end - region.start) + 1;
    }
}

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

int Tx::get_cds_start() {
    char fwd = '+';
    if (get_strand() == fwd) {
        return cds_min;
    } else {
        return cds_max;
    }
}

int Tx::get_cds_end() {
    char fwd = '+';
    if (get_strand() == fwd) {
        return cds_max;
    } else {
        return cds_min;
    }
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

// include a sort operator for Region structs. This allows quick searching
// of vectors of Regions
static bool compareRegion(const Region& a, int b) {
    return a.start < b;
}

// find the index of the closest exon to a chromosome position
int Tx::closest_exon_num(int pos) {
    return closest_exon_num(pos, exons);
}

// find the index for the closest exon/CDS to a chromosome position
//
// @param position chromosomal nucleotide position
// @returns number of exon containing the position
int Tx::closest_exon_num(int pos, std::vector<Region> & group) {
    size_t idx = std::lower_bound(group.begin(), group.end(), pos, compareRegion) - group.begin();
    size_t size = group.size();
    if (idx == 0) {
        return idx;
    } else if ((idx == size) && (pos > group[idx - 1].start)) {
        return idx - 1;
    }
    
    // get the exon, and the preceding exon (lower bound finds exon after pos)
    Region a = group[idx - 1];
    Region b = group[idx];
    
    // if the position is within an exon, return that exon number
    if ((pos >= a.start) && (pos <= a.end)) {
        return idx - 1;
    } else if ((pos >= b.start) && (pos <= b.end)) {
        return idx;
    }
    
    // pos lies after the start of a. a is upstream, b is downstream, find closest.
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
    if (cds.size() == 0) {
        return false;
    }
    int idx = closest_exon_num(pos, cds);
    Region region = cds[idx];
    return (pos >= region.start) && (pos <= region.end);
}

// shift an intronic site to the nearest exon boundary (but track the offset)
CDS_coords Tx::to_closest_exon(int pos, Region exon) {
    // Region exon = get_closest_exon(pos);
    int site;
    if (std::abs(exon.start - pos) < std::abs(exon.end - pos)) {
        site = exon.start;
    } else {
        site = exon.end;
    }
    
    bool fwd = get_strand() == '+';
    int offset = (fwd) ? pos - site : site - pos;
    return CDS_coords {site, offset};
}

// find the distance in coding bp to the CDS start
//
// @param pos first chromosome position
// @returns CDS_coords, which contains distance in coding sequence base
//    pairs between two positions and an offset for distance into the
//    intron
CDS_coords Tx::get_coding_distance(int pos) {
    int cds_start = get_cds_start();
    
    int first = closest_exon_num(cds_start);
    int alt = closest_exon_num(pos);

    int offset = 0;
    if ( !is_exonic(pos, exons[alt]) ) {
        CDS_coords site = to_closest_exon(pos, exons[alt]);
        offset = site.offset;
        pos = site.position;
    }
    
    bool fwd = get_strand() == '+';
    bool after = alt > first;
    bool positive = (fwd & after) | (!fwd & !after);
    // use positive coding position if the transcript is either:
    //   a) on the + strand and the site is 3' of CDS start
    //   b) on the - strand and the site is 5' of CDS start
    int delta = 0;
    if (first == alt) {
        delta = (fwd) ? pos - cds_start : cds_start - pos ;
    } else {
        int dist_1;
        int dist_2;
        if (after) {
            dist_1 = (fwd) ? exons[first].end - cds_start : cds_start - exons[first].end;
            dist_2 = (fwd) ? (pos - exons[alt].start) + 1 : (exons[alt].start - pos) - 1;
        } else {
            dist_1 = (fwd) ? exons[first].start - cds_start : cds_start - exons[first].start;
            dist_2 = (fwd) ? (pos - exons[alt].end) - 1 : (exons[alt].end - pos) + 1;
        }
        
        int lo = std::min(first, alt) + 1;
        int hi = std::max(first, alt);
        delta = dist_1 + dist_2;
        for (int i=lo; i < hi; i++) {
            int left = exons[i].start;
            int right = exons[i].end;
            delta += (positive) ? (right - left) + 1 : (left - right) - 1;
        }
    }
    
    return CDS_coords {delta, offset};
}

// cache the exon boundary positions as CDS coordinates. 
//
// Rather than recalculate the CDS positions each time we want to convert a 
// chromosome position to a CDS position, we calculate them once and cache 
// the coordinates in a dictionary for rapid access.
void Tx::_cache_exon_cds_positions() {
    exon_to_cds.clear();
    
    for (auto &region : cds) {
        // get exon boundaries as CDS distance from the start site
        CDS_coords start = get_coding_distance(region.start);
        CDS_coords end = get_coding_distance(region.end);
        
        // cache the CDS positions of the exon boundaries
        exon_to_cds[region.start] = std::abs(start.position);
        exon_to_cds[region.end] = std::abs(end.position);
    }
}

// figure out the chromosome position of a CDS site
//
// @param cds_position position of a variant in CDS.
// @returns chromosome bp position of the CDS site
int Tx::get_position_on_chrom(int cds_position, int offset) {
    // cache the exon boundaries in CDS distances
    if (exon_to_cds.empty()) { _cache_exon_cds_positions(); }
    
    int start;
    int end;
    int start_cds;
    int end_cds;
    
    // quickly find the exon containing the CDS position
    for (auto &region : cds) {
        start = region.start;
        end = region.end;
        
        start_cds = std::min(exon_to_cds[start], exon_to_cds[end]);
        end_cds = std::max(exon_to_cds[start], exon_to_cds[end]);
        
        if (start_cds <= cds_position && cds_position <= end_cds) {
            // convert the CDS position to a chromosomal coordinate
            char fwd = '+';
            if (get_strand() == fwd) {
                return start + (cds_position - start_cds) + offset;
            } else {
                return start + (end_cds - cds_position) - offset;
            }
        }
    }
    
    throw std::invalid_argument( "position not in CDS" );
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
    cds_sequence = std::move(cds_dna);
}

// add and process the genomic sequence into a CDS sequence
//
// @param gdna string for genomic sequence. If the transcript strand is '-',
//    then the DNA sequence will be for the - strand, so we need to
//    reorient the DNA to the + strand.
void Tx::add_genomic_sequence(std::string gdna) {
    char fwd = '+';
    if (get_strand() != fwd) {
        // orient the DNA sequence to the + strand.
        gdna = reverse_complement(gdna);
    }
    
    genomic_sequence = std::move(gdna);
    
    std::string cds_seq;
    for (auto &region : cds) {
        int x = std::abs(region.start - get_start()) + gdna_offset;
        int len = (region.end - region.start) + 1;
        cds_seq += genomic_sequence.substr(x, len);
    }
    if (get_strand() != fwd) {
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
    
    _fix_cds_length();
    
    cds_length = 0;
    for (auto region : cds) {
        cds_length += (region.end - region.start) + 1;
    }
}

// correct the coding sequence of a transcript, if it misses bases.
//
// Some transcripts don't cover the full length of a gene, they terminate
// in the middle of an amino acid. This is rare, and it is rarer that we
// select these transcripts for analysis. TNS3 is an example, where the
// de novos in the gene can only be contained within a transcript that is
// incomplete at the 3' end. This causes problems when we try to translate
// the final codon.
//
// We simply extend the coding sequence 1-2 bases to make a complete codon.
void Tx::_fix_cds_length() {
    if (cds.empty()) {
        return;
    }
    exon_to_cds.clear();
    int diff = cds_sequence.size() % 3;
    int end = get_cds_end();

    size_t last = cds.size() - 1;

    // genomic_sequence is stored in + strand orientation, so the extra bases
    // come from just past the 3' CDS end on the + strand, or just before it
    // (then reverse complemented) on the - strand
    if (diff != 0) {
        diff = 3 - diff;
        
        bool fwd = get_strand() == '+';
        int start_bp = (fwd ? end + 1 : end - diff) - get_start() + gdna_offset;
        if (start_bp < 0 || start_bp + diff > (int) genomic_sequence.size()) {
            // the genomic sequence doesn't extend far enough to complete the codon
            return;
        }
        
        if (fwd) {
            cds[last] = Region {cds[last].start, cds[last].end + diff};
            cds_sequence += genomic_sequence.substr(start_bp, diff);
        } else {
            cds[0] = Region {cds[0].start - diff, cds[0].end};
            cds_sequence += reverse_complement(genomic_sequence.substr(start_bp, diff));
        }
    }
    
    cds_min = cds[0].start;
    cds_max = cds[last].end;
    
    // shifting the CDS coordinates can, very infrequently, shift the CDS
    // beyond the exon coordinates. This is only a problem if we take the union
    // transcript coordinates, but should be accounted for.
    if (cds_min < exons[0].start) {
        exons[0] = Region {exons[0].start - diff, exons[0].end};
    }
    
    last = exons.size() - 1;
    if (cds_max > exons[last].end) {
        exons[last] = Region {exons[last].start, exons[last].end + diff};
    }
}

// lookup table of complementary bases, including IUPAC ambiguity codes. Any
// other character complements to 'N'.
static std::array<char, 256> make_complements() {
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
static char complement(char base) {
    static const std::array<char, 256> complements = make_complements();
    return complements[(uint8_t) base];
}

// reverse complement a DNA or RNA sequence
std::string Tx::reverse_complement(const std::string & seq) {
    std::string revcomp(seq.rbegin(), seq.rend());
    for (auto &base : revcomp) {
        base = complement(base);
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
    int end = start + length;
    return get_seq_in_region(start, end);
}

std::string Tx::get_seq_in_region(int start, int end) {
    if (start < (get_start() - gdna_offset) || end > (get_end() + gdna_offset)) {
        throw std::invalid_argument( "sequence position not in gene range" );
    }
    int pos = start - get_start() + gdna_offset;
    int length = end - start;
    std::string seq = genomic_sequence.substr(pos, length);
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    return seq;
}

// get the codon sequence for a given codon_number
std::string Tx::get_codon_sequence(int codon) { 
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
    if (aa_code.count(codon) == 0) {
        std::string msg = "cannot translate codon: " + codon;
        throw std::invalid_argument(msg);
    }
    return aa_code[codon];
}

// translate a DNA codon to a single character amino acid
//
// @param seq codon sequence, or longer DNA sequences.
// @returns translated amino acid sequence as single character codes.
std::string Tx::translate(std::string seq) {
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    
    std::string protein;
    protein.resize(seq.size() / 3);
    
    for (std::uint32_t i=0, j=0; i<seq.size(); i=i+3, j++) {
        protein[j] = translate_codon(seq.substr(i, 3));
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
    
    if ((site.position < 0) | (site.position >= cds_length)) {
        throw std::invalid_argument( "position not inside CDS region" );
    }
    
    // define the default values for the codon positions. We check these later
    // in python and convert to None if they are still the defaults.
    int codon_number = -9999999;
    int intra_codon = -1;
    std::string codon_seq = "";
    char initial_aa='0';
    
    if (in_coding) {
        codon_number = get_codon_number_for_cds_position(site.position);
        intra_codon = get_position_within_codon(site.position);
        codon_seq = get_codon_sequence(codon_number);
        initial_aa = translate_codon(codon_seq);
    }
    
    return Codon {site.position, codon_seq, intra_codon, codon_number,
        initial_aa, site.offset};
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

// check for overlap of two ranges
static bool overlap(const int & start1, const int & end1, const int & start2, const int & end2) {
    return (end1 >= start2) && (end2 >= start1);
}

// find the index of a region overlapping a genomic range, checking the region
// closest to the range start, then its neighbours. Returns -1 if none overlap.
int Tx::overlapping_region_num(int start, int end, std::vector<Region> & group) {
    if (group.empty()) {
        return -1;
    }
    int idx = closest_exon_num(start, group);
    int prev = std::max(0, idx - 1);
    int post = std::min((int) group.size() - 1, idx + 1);
    for (int i : {idx, prev, post}) {
        if (overlap(start, end, group[i].start, group[i].end)) {
            return i;
        }
    }
    return -1;
}

// checks if a genomic range overlaps any CDS range
bool Tx::overlaps_cds(int start, int end) {
    return overlapping_region_num(start, end, cds) >= 0;
}

// checks if a genomic range overlaps the 3 exonic bases next to an intron.
// Introns under 10 bp are treated as exonic, so don't have a splice region.
bool Tx::in_exonic_splice_region(int start, int end) {
    int idx = overlapping_region_num(start, end, exons);
    if (idx < 0) {
        return false;
    }
    Region exon = exons[idx];
    if ((idx > 0) && (exon.start - exons[idx - 1].end >= 10)
            && overlap(start, end, exon.start, exon.start + 2)) {
        return true;
    }
    if ((idx < (int) exons.size() - 1) && (exons[idx + 1].start - exon.end >= 10)
            && overlap(start, end, exon.end - 2, exon.end)) {
        return true;
    }
    return false;
}

// consequences in order of severity, following Ensembl VEP
static const std::vector<std::string> cq_order = {"splice_acceptor_variant",
    "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost",
    "start_lost", "inframe_insertion", "inframe_deletion", "missense_variant",
    "protein_altering_variant", "splice_donor_5th_base_variant",
    "splice_region_variant", "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant", "start_retained_variant",
    "stop_retained_variant", "synonymous_variant", "coding_sequence_variant",
    "5_prime_UTR_variant", "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant", "intron_variant",
    "upstream_gene_variant", "downstream_gene_variant", "intergenic_variant"};

static size_t severity(const std::string & cq) {
    return std::find(cq_order.begin(), cq_order.end(), cq) - cq_order.begin();
}

static std::string most_severe(const std::string & a, const std::string & b) {
    return (severity(a) <= severity(b)) ? a : b;
}

// consequence for intronic sites from near to far bases from the splice donor
// (at the 5' end of the intron)
static std::string donor_cq(int near, int far) {
    if (overlap(near, far, 1, 2)) { return "splice_donor_variant"; }
    if (overlap(near, far, 5, 5)) { return "splice_donor_5th_base_variant"; }
    if (overlap(near, far, 7, 8)) { return "splice_region_variant"; }
    if (overlap(near, far, 3, 6)) { return "splice_donor_region_variant"; }
    return "intron_variant";
}

// consequence for intronic sites from near to far bases from the splice
// acceptor (at the 3' end of the intron)
static std::string acceptor_cq(int near, int far) {
    if (overlap(near, far, 1, 2)) { return "splice_acceptor_variant"; }
    if (overlap(near, far, 3, 8)) { return "splice_region_variant"; }
    if (overlap(near, far, 3, 17)) { return "splice_polypyrimidine_tract_variant"; }
    return "intron_variant";
}

// consequence for an intronic insertion, where near is the distance from the
// splice donor of the closer base flanking the insertion
static std::string donor_insertion_cq(int near) {
    if (near == 1) { return "splice_donor_variant"; }
    if (near >= 3 && near <= 5) { return "splice_donor_region_variant"; }
    if (near == 2 || near == 6 || near == 7) { return "splice_region_variant"; }
    return "intron_variant";
}

// consequence for an intronic insertion, where near is the distance from the
// splice acceptor of the closer base flanking the insertion
static std::string acceptor_insertion_cq(int near) {
    if (near == 1) { return "splice_acceptor_variant"; }
    if (near >= 2 && near <= 7) { return "splice_region_variant"; }
    if (near >= 8 && near <= 16) { return "splice_polypyrimidine_tract_variant"; }
    return "intron_variant";
}

// consequence of replacing a reference peptide with an altered peptide, if
// that gains, loses or retains a stop, or loses the start. Otherwise returns
// an empty string. As in VEP, the start is only lost if the altered peptide
// neither starts nor ends with the reference peptide, so deleting whole codons
// including the start codon is not start_lost.
static std::string peptide_cq(const std::string & ref_aa, const std::string & alt_aa, bool at_start) {
    size_t ref_stop = ref_aa.find('*');
    size_t alt_stop = alt_aa.find('*');
    if (ref_stop == std::string::npos && alt_stop != std::string::npos) {
        return "stop_gained";
    } else if (ref_stop != std::string::npos && alt_stop == std::string::npos) {
        return "stop_lost";
    } else if (ref_stop != std::string::npos
            && ref_aa.compare(0, ref_stop + 1, alt_aa, 0, alt_stop + 1) == 0) {
        return "stop_retained_variant";
    } else if (at_start && !alt_aa.empty()) {
        bool starts = alt_aa.compare(0, ref_aa.size(), ref_aa) == 0;
        bool ends = alt_aa.size() >= ref_aa.size()
            && alt_aa.compare(alt_aa.size() - ref_aa.size(), ref_aa.size(), ref_aa) == 0;
        if (!starts && !ends) {
            return "start_lost";
        }
    }
    return "";
}

static bool is_stop_change(const std::string & cq) {
    return cq == "stop_gained" || cq == "stop_lost" || cq == "stop_retained_variant";
}

std::string Tx::outside_gene_cq(int start, int end) {
    bool fwd = get_strand() == '+';
    if ((end < tx_start) & ((tx_start - end) <= 5000)) {
        if (fwd) {
            return "upstream_gene_variant";
        } else {
            return "downstream_gene_variant";
        }
    } else if ((start > tx_end) & ((start - tx_end) <= 5000)) {
        if (fwd) {
            return "downstream_gene_variant";
        } else {
            return "upstream_gene_variant";
        }
    }
    return "intergenic_variant";
}

// consequence for a site in an intron under 10 bp, or an empty string for
// sites elsewhere. Introns this short are likely annotation artefacts, so the
// site is treated as exonic.
std::string Tx::short_intron_cq(int start, int end) {
    auto idx = closest_exon_num(start);
    auto exon = exons[idx];
    int intron_length = 100;
    if ((exon.end < start) & (idx < ((int) exons.size() - 1))) {
        auto other = exons[idx + 1];
        intron_length = other.start - exon.end;
    } else if ((end < exon.start) & (idx > 0)) {
        auto other = exons[idx - 1];
        intron_length = exon.start - other.end;
    }
    if (intron_length >= 10) {
        return "";
    }
    bool fwd = get_strand() == '+';
    if (cds.empty()) {
        return "non_coding_transcript_exon_variant";
    } else if (end < cds_min) {
        return fwd ? "5_prime_UTR_variant" : "3_prime_UTR_variant";
    } else if (start > cds_max) {
        return fwd ? "3_prime_UTR_variant" : "5_prime_UTR_variant";
    }
    return "coding_sequence_variant";
}

// consequence for sites within an intron, checking both the splice donor (at
// the 5' end of the intron) and the splice acceptor (at the 3' end). For
// insertions, start and end are the bases either side of the inserted bases.
std::string Tx::intron_cq(int start, int end, bool insertion) {
    int idx = closest_exon_num(start);
    if (exons[idx].end >= start) {
        // the site precedes the closest exon
        idx -= 1;
    }
    if (idx < 0 || idx >= (int) exons.size() - 1) {
        return "intron_variant";
    }
    
    // distances into the intron, from the exon before, and the exon after
    int lower = exons[idx].end;
    int upper = exons[idx + 1].start;
    bool fwd = get_strand() == '+';
    int donor_near = fwd ? start - lower : upper - end;
    int donor_far = fwd ? end - lower : upper - start;
    int acceptor_near = fwd ? upper - end : start - lower;
    int acceptor_far = fwd ? upper - start : end - lower;
    if (insertion) {
        return most_severe(donor_insertion_cq(donor_near), acceptor_insertion_cq(acceptor_near));
    }
    return most_severe(donor_cq(donor_near, donor_far), acceptor_cq(acceptor_near, acceptor_far));
}

// figure out consequence for a site not in coding region
std::string Tx::intronic_cq(int start, int end) {
    std::string short_intron = short_intron_cq(start, end);
    if (!short_intron.empty()) {
        return short_intron;
    }
    
    bool fwd = get_strand() == '+';
    // TODO: figure out for variants spanning an exon
    CDS_coords start_coord = get_coding_distance(start);
    CDS_coords end_coord = get_coding_distance(end);
    
    // variants running from an exon into an intron disrupt that splice site
    if ((start_coord.offset == 0) & (end_coord.offset != 0) & (end < exons.back().end)) {
        return fwd ? "splice_donor_variant" : "splice_acceptor_variant";
    } else if ((start_coord.offset != 0) & (end_coord.offset == 0) & (start > exons[0].start)) {
        return fwd ? "splice_acceptor_variant" : "splice_donor_variant";
    }
    
    bool intronic = (start_coord.offset != 0) & (end_coord.offset != 0);
    if (!intronic) {
        std::string cq;
        if (cds.empty()) {
            cq = "non_coding_transcript_exon_variant";
        } else if (end_coord.position <= 0) {
            cq = "5_prime_UTR_variant";
        } else if (start_coord.position >= cds_length) {
            cq = "3_prime_UTR_variant";
        }
        if (!cq.empty()) {
            return in_exonic_splice_region(start, end) ? "splice_region_variant" : cq;
        }
    }
    
    // TODO: indels near the splice region can duplicate the local sequence.
    // TODO: In those cases, the splice site might remain unaltered. This means
    // TODO: the variant location might need to be adjusted. For this, we'd need 
    // TODO: to check the ref genome vs the replacement, to adjust indel coords.
    return intron_cq(start, end, false);
}

// up to n bases of the 3' UTR following the CDS, in transcript orientation
std::string Tx::utr3_sequence(int n) {
    std::string seq;
    bool fwd = get_strand() == '+';
    int step = fwd ? 1 : -1;
    int cds_end = get_cds_end();
    int first = closest_exon_num(cds_end);
    for (int i = first; i >= 0 && i < (int) exons.size() && (int) seq.size() < n; i += step) {
        int from = (i == first) ? cds_end + step : (fwd ? exons[i].start : exons[i].end);
        int to = fwd ? exons[i].end : exons[i].start;
        for (int pos = from; (to - pos) * step >= 0 && (int) seq.size() < n; pos += step) {
            int x = pos - get_start() + gdna_offset;
            if (x < 0 || x >= (int) genomic_sequence.size()) {
                return seq;
            }
            seq += fwd ? genomic_sequence[x] : complement(genomic_sequence[x]);
        }
    }
    return seq;
}

// translate the codons covering CDS positions lo to hi, before and after
// replacing those bases with alt_t (in transcript orientation). For
// insertions, lo is the position the bases are inserted before, and hi is
// lo - 1. After a frameshift in the stop codon, the last alt codon is
// completed with the following bases, to check if the stop is retained, but
// otherwise incomplete codons aren't translated. Codons running past the CDS
// end use the 3' UTR sequence.
//
// @returns false if the sequence isn't available
bool Tx::translate_change(int lo, int hi, const std::string & alt_t,
        std::string & ref_aa, std::string & alt_aa) {
    if (lo < 0 || cds_sequence.empty() || (int) cds_sequence.size() != cds_length) {
        return false;
    }
    int last = std::max(lo, hi);
    int first = lo - lo % 3;
    int length = last - last % 3 + 3 - first;
    if (hi < lo && lo % 3 == 0 && last < cds_length - 3) {
        // insertions between codons leave the neighbouring codons intact
        length = 0;
    }
    int alt_length = length - (hi - lo + 1) + (int) alt_t.size();
    int needed = first + length;
    if (last >= cds_length - 3) {
        needed += (3 - alt_length % 3) % 3;
    }
    std::string seq = cds_sequence;
    if (needed > cds_length) {
        seq += utr3_sequence(needed - cds_length);
    }
    if (first + length > (int) seq.size()) {
        return false;
    }
    std::string ref_seq = seq.substr(first, length);
    std::string alt_seq = ref_seq;
    alt_seq.replace(lo - first, hi - lo + 1, alt_t);
    alt_seq += seq.substr(first + length, needed - first - length);
    alt_seq.resize(alt_seq.size() - alt_seq.size() % 3);
    try {
        ref_aa = translate(ref_seq);
        alt_aa = translate(alt_seq);
    } catch (const std::invalid_argument &) {
        return false;
    }
    return true;
}

// consequence for bases inserted before CDS position k
std::string Tx::coding_insertion_cq(int k, const std::string & ins_t) {
    std::string ref_aa;
    std::string alt_aa;
    std::string change;
    if (translate_change(k, k - 1, ins_t, ref_aa, alt_aa)) {
        change = peptide_cq(ref_aa, alt_aa, k < 3);
    }
    if (ins_t.size() % 3 != 0) {
        return is_stop_change(change) ? change : "frameshift_variant";
    } else if (!change.empty()) {
        return change;
    }
    // insertions within a codon can change the amino acid at that codon
    bool retained = (alt_aa.compare(0, ref_aa.size(), ref_aa) == 0)
        || (alt_aa.size() >= ref_aa.size()
            && alt_aa.compare(alt_aa.size() - ref_aa.size(), ref_aa.size(), ref_aa) == 0);
    return retained ? "inframe_insertion" : "protein_altering_variant";
}

// consequence for bases inserted between the left base and the next base
//
// @param left position of the base before the inserted bases
// @param ins inserted bases, on the + strand
std::string Tx::insertion_cq(int left, const std::string & ins) {
    int right = left + 1;
    if (right <= tx_start) {
        return outside_gene_cq(left, left);
    } else if (left >= tx_end) {
        return outside_gene_cq(right, right);
    }
    
    bool left_exonic = is_exonic(left);
    bool right_exonic = is_exonic(right);
    if (!left_exonic && !right_exonic) {
        std::string short_intron = short_intron_cq(left, right);
        return short_intron.empty() ? intron_cq(left, right, true) : short_intron;
    }
    
    // insertions at an exon/intron junction count as being in the exon
    bool fwd = get_strand() == '+';
    std::string cq;
    if (cds.empty()) {
        cq = "non_coding_transcript_exon_variant";
    } else {
        // find the CDS position the bases are inserted before
        int before = fwd ? left : right;
        int after = fwd ? right : left;
        bool before_exonic = fwd ? left_exonic : right_exonic;
        int k = before_exonic ? get_coding_distance(before).position + 1
            : get_coding_distance(after).position;
        if (k <= 0) {
            cq = "5_prime_UTR_variant";
        } else if (k >= cds_length) {
            cq = "3_prime_UTR_variant";
        } else {
            cq = coding_insertion_cq(k, fwd ? ins : reverse_complement(ins));
        }
    }
    
    // insertions within the 3 exonic bases next to an intron, or at the
    // exon/intron junction, are in the splice region
    bool region = (!left_exonic || in_exonic_splice_region(left, left))
        && (!right_exonic || in_exonic_splice_region(right, right));
    return region ? most_severe(cq, "splice_region_variant") : cq;
}

// check if an indel within a single exon gains, loses or retains a stop, or
// loses the start, by translating the affected codons before and after.
//
// @returns consequence, or an empty string if none of those apply, or the
//    sequence isn't available
std::string Tx::indel_codon_cq(int start, int end, const std::string & ref, const std::string & alt) {
    CDS_coords start_coord = get_coding_distance(start);
    CDS_coords end_coord = get_coding_distance(end);
    int lo = std::min(start_coord.position, end_coord.position);
    int hi = std::max(start_coord.position, end_coord.position);
    if (start_coord.offset != 0 || end_coord.offset != 0 || (hi - lo + 1) != (int) ref.size()) {
        return "";
    }
    if (lo < 0) {
        return "start_lost";
    }
    
    std::string ref_aa;
    std::string alt_aa;
    std::string alt_t = (get_strand() == '+') ? alt : reverse_complement(alt);
    if (!translate_change(lo, hi, alt_t, ref_aa, alt_aa)) {
        return "";
    }
    std::string cq = peptide_cq(ref_aa, alt_aa, lo < 3);
    if (cq.empty() && lo < 3) {
        // as in VEP, the start is also lost if the CDS no longer begins with ATG
        std::string altered = cds_sequence.substr(0, hi + 4);
        altered.replace(lo, hi - lo + 1, alt_t);
        if (altered.compare(0, 3, "ATG") != 0) {
            return "start_lost";
        }
    }
    return cq;
}

// check consequences for an indel with some CDS overlap
std::string Tx::indel_cq(int start, int end, std::string ref, std::string alt) {
    // find the overlapped exon, then check how much of the indel lies in that
    int idx = overlapping_region_num(start, end, exons);
    if (idx < 0) {
        throw (std::invalid_argument("cannot find matching exon for indel"));
    }
    Region exon = exons[idx];
    
    // deletions running into an intron disrupt the splice site, which is more
    // severe than any frameshift or inframe change. Deletions running past
    // the transcript ends are trimmed to the exonic part instead.
    bool fwd = get_strand() == '+';
    bool deletion = ref.size() > alt.size();
    if (deletion & (start < exon.start)) {
        if (idx > 0) {
            return fwd ? "splice_acceptor_variant" : "splice_donor_variant";
        }
        ref = ref.substr(exon.start - start);
        alt = "";
        start = exon.start;
    } else if (deletion & (end > exon.end)) {
        if (idx < (int) exons.size() - 1) {
            return fwd ? "splice_donor_variant" : "splice_acceptor_variant";
        }
        ref = ref.substr(0, (exon.end - start) + 1);
        alt = "";
        end = exon.end;
    }

    int ref_len = (int) ref.size();
    int alt_len = (int) alt.size();
    bool inframe = (std::abs(ref_len - alt_len) % 3) == 0;
    std::string cq = indel_codon_cq(start, end, ref, alt);
    bool into_utr5 = fwd ? start < get_cds_start() : end > get_cds_start();
    if (!inframe && into_utr5) {
        // frameshifts must lie within the CDS, but this disrupts the start
        cq = "start_lost";
    } else if (!inframe) {
        cq = is_stop_change(cq) ? cq : "frameshift_variant";
    } else if (cq.empty()) {
        // TODO: account for altered splice boundaries
        cq = deletion ? "inframe_deletion" : "inframe_insertion";
    }
    return in_exonic_splice_region(start, end) ? most_severe(cq, "splice_region_variant") : cq;
}

// consequence of changing a single codon
static std::string codon_cq(char initial_aa, char mutated_aa, int codon_number, int last_codon) {
    if (initial_aa != '*' && mutated_aa == '*') { return "stop_gained"; }
    if (initial_aa == '*' && mutated_aa != '*') { return "stop_lost"; }
    
    if (initial_aa != mutated_aa) {
        if (codon_number == 0) {
            return "start_lost";
        }
        return "missense_variant";
    }
    
    if (codon_number == 0) {
        return "start_retained_variant";
    } else if (codon_number == last_codon && initial_aa == '*') {
        return "stop_retained_variant";
    }
    return "synonymous_variant";
}

// consequence of a SNV or MNV which overlaps the CDS
//
// Each substituted base within the CDS is applied to its codon, and the most
// severe consequence across the altered codons is returned. Bases outside the
// CDS are skipped.
std::string Tx::coding_cq(int start, int end, const std::string & alt) {
    if (cds_sequence.empty()) {
        throw std::invalid_argument(get_name() + " lacks genomic sequence (e.g. "
            "from a genome fasta), which is needed for coding consequences");
    }
    bool fwd = get_strand() == '+';
    std::map<int, Codon> codons;
    for (int pos = start; pos <= end; pos++) {
        if (!in_coding_region(pos)) {
            continue;
        }
        Codon codon = get_codon_info(pos);
        auto it = codons.emplace(codon.codon_number, codon).first;
        char base = alt[pos - start];
        if (!fwd) {
            base = complement(base);
        }
        it->second.codon_seq[codon.intra_codon] = base;
    }
    
    if (codons.empty()) {
        return "coding_sequence_variant";
    }
    
    int last_codon = get_codon_number_for_cds_position(cds_length - 1);
    std::string worst = "synonymous_variant";
    for (auto &item : codons) {
        Codon & codon = item.second;
        char mutated_aa = translate_codon(codon.codon_seq);
        worst = most_severe(worst, codon_cq(codon.initial_aa, mutated_aa, codon.codon_number, last_codon));
    }
    if (in_exonic_splice_region(start, end)) {
        worst = most_severe(worst, "splice_region_variant");
    }
    return worst;
}

int min_len(const std::string & a, const std::string & b) {
    return (int)std::min(a.size(), b.size());
}

// realign alt against ref sequence, and strip down to alternate bases
void Tx::trim_alleles(int& start, int&end, std::string& ref, std::string& alt) {
    // trim the front of the alleles
    while ((min_len(ref, alt) > 0) && (ref.size() > 1) && (alt[0] == ref[0])) {
        ref.erase(0, 1);
        alt.erase(0, 1);
        start += 1;
    }
    
    // trim the back of the alleles
    while ((min_len(ref, alt) > 0) && (ref.size() > 1) && (alt.back() == ref.back())) {
        ref.pop_back();
        alt.pop_back();
        end -= 1;
    }
}

std::string Tx::consequence(int pos, std::string ref, std::string alt) {
    int start = pos;
    int end = pos + (int)ref.size() - 1;
    if (ref.empty()) {
        throw std::invalid_argument("reference allele is empty");
    }
    
    // TODO: need to realign alt against ref sequence, and strip down to alternate bases
    trim_alleles(start, end, ref, alt);
    
    // insertions lie between the anchor base and its neighbour
    if (ref.size() == 1 && alt.size() > 1 && (alt[0] == ref[0] || alt.back() == ref[0])) {
        if (alt[0] == ref[0]) {
            return insertion_cq(start, alt.substr(1));
        }
        return insertion_cq(start - 1, alt.substr(0, alt.size() - 1));
    }
    
    if ((end < tx_start) | (start > tx_end)) {
        return outside_gene_cq(start, end);
    }
    
    if (!overlaps_cds(start, end)) {
        return intronic_cq(start, end);
    }
    
    // at this point we know that the site lies within the gene boundary, and
    // at least partially overlaps the CDS
    if (ref.size() != alt.size()) {
        return indel_cq(start, end, ref, alt);
    }
    
    return coding_cq(start, end, alt);
}

} // namespace