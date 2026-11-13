#include <algorithm>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

#include "tx.h"
#include "bases.h"

// predicting the consequences of variants on a transcript, following Ensembl VEP

namespace gencode {

// introns are treated as exonic (as likely annotation artefacts) if the next
// exon starts fewer than this many bases after the previous exon ends
static const int min_intron_gap = 10;

// sites within this distance of a transcript are upstream or downstream of it
static const int flank_distance = 5000;

// checks if a genomic range overlaps the 3 exonic bases next to an intron.
// Short introns are treated as exonic, so don't have a splice region.
bool Tx::in_exonic_splice_region(int start, int end) {
    int idx = overlapping_region_num(start, end, exons);
    if (idx < 0) {
        return false;
    }
    Region exon = exons[idx];
    if ((idx > 0) && (exon.start - exons[idx - 1].end >= min_intron_gap)
            && overlap(start, end, exon.start, exon.start + 2)) {
        return true;
    }
    if ((idx < (int) exons.size() - 1) && (exons[idx + 1].start - exon.end >= min_intron_gap)
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

// consequence for a genomic range running out of an overlapped exon into an
// adjacent intron, which disrupts that splice site. Returns an empty string if
// the range stays within the exon, or only runs past the transcript ends.
//
// @param idx index of an exon overlapping the range
std::string Tx::splice_boundary_cq(int start, int end, int idx) {
    bool fwd = is_forward();
    std::string cq;
    if (idx > 0 && start < exons[idx].start) {
        cq = fwd ? "splice_acceptor_variant" : "splice_donor_variant";
    }
    if (idx >= 0 && idx < (int) exons.size() - 1 && end > exons[idx].end) {
        cq = most_severe(cq, fwd ? "splice_donor_variant" : "splice_acceptor_variant");
    }
    return cq;
}

// consequence for sites outside the transcript, which are upstream or
// downstream if nearby
std::string Tx::outside_gene_cq(int start, int end) {
    bool before = end < tx_start;
    int distance = before ? tx_start - end : start - tx_end;
    if (distance > flank_distance) {
        return "intergenic_variant";
    }
    return (before == is_forward()) ? "upstream_gene_variant" : "downstream_gene_variant";
}

// consequence for a site in a short intron (see min_intron_gap), or an empty
// string for sites elsewhere. The site is treated as exonic.
std::string Tx::short_intron_cq(int start, int end) {
    auto idx = closest_exon_num(start);
    auto exon = exons[idx];
    int gap = min_intron_gap;
    if ((exon.end < start) && (idx < ((int) exons.size() - 1))) {
        gap = exons[idx + 1].start - exon.end;
    } else if ((end < exon.start) && (idx > 0)) {
        gap = exon.start - exons[idx - 1].end;
    }
    if (gap >= min_intron_gap) {
        return "";
    }
    bool fwd = is_forward();
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
    bool fwd = is_forward();
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
    
    // variants running from an exon into an intron disrupt that splice site
    std::string splice = splice_boundary_cq(start, end, overlapping_region_num(start, end, exons));
    if (!splice.empty()) {
        return splice;
    }
    
    CDS_coords start_coord = get_coding_distance(start);
    CDS_coords end_coord = get_coding_distance(end);
    
    bool intronic = (start_coord.offset != 0) && (end_coord.offset != 0);
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
    bool fwd = is_forward();
    int step = fwd ? 1 : -1;
    int cds_end = get_cds_end();
    int first = closest_exon_num(cds_end);
    for (int i = first; i >= 0 && i < (int) exons.size() && (int) seq.size() < n; i += step) {
        int from = (i == first) ? cds_end + step : (fwd ? exons[i].start : exons[i].end);
        int to = fwd ? exons[i].end : exons[i].start;
        for (int pos = from; (to - pos) * step >= 0 && (int) seq.size() < n; pos += step) {
            int x = genomic_index(pos);
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
    bool fwd = is_forward();
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
    std::string alt_t = is_forward() ? alt : reverse_complement(alt);
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
    
    // indels running into an intron disrupt the splice site, which is more
    // severe than any frameshift or inframe change. Deletions running past
    // the transcript ends are trimmed to the exonic part instead.
    bool fwd = is_forward();
    bool deletion = ref.size() > alt.size();
    std::string splice = splice_boundary_cq(start, end, idx);
    if (!splice.empty()) {
        return splice;
    }
    if (alt.empty() && start < exon.start) {
        ref = ref.substr(exon.start - start);
        start = exon.start;
    }
    if (alt.empty() && end > exon.end) {
        ref = ref.substr(0, (exon.end - start) + 1);
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
// CDS, or in an incomplete final codon, are skipped.
std::string Tx::coding_cq(int start, int end, const std::string & alt) {
    if (cds_sequence.empty()) {
        throw std::invalid_argument(get_name() + " lacks genomic sequence (e.g. "
            "from a genome fasta), which is needed for coding consequences");
    }
    bool fwd = is_forward();
    int partial = cds_sequence.size() % 3;
    int complete = (int) cds_sequence.size() - partial;
    // CDS positions change monotonically with chromosome position, so bases
    // in the same codon are consecutive
    std::vector<Codon> codons;
    for (int pos = start; pos <= end; pos++) {
        if (!in_coding_region(pos)) {
            continue;
        }
        if (partial != 0 && get_coding_distance(pos).position >= complete) {
            continue;
        }
        Codon codon = get_codon_info(pos);
        if (codons.empty() || codons.back().codon_number != codon.codon_number) {
            codons.push_back(codon);
        }
        char base = alt[pos - start];
        if (!fwd) {
            base = complement(base);
        }
        codons.back().codon_seq[codon.intra_codon] = base;
    }
    
    if (codons.empty()) {
        return "coding_sequence_variant";
    }
    
    int last_codon = get_codon_number_for_cds_position(cds_length - 1);
    std::string worst = "synonymous_variant";
    for (auto &codon : codons) {
        char mutated_aa = translate_codon(codon.codon_seq);
        worst = most_severe(worst, codon_cq(codon.initial_aa, mutated_aa, codon.codon_number, last_codon));
    }
    if (in_exonic_splice_region(start, end)) {
        worst = most_severe(worst, "splice_region_variant");
    }
    
    // changes running from the exon into an intron disrupt the splice site.
    // Single base changes in the CDS can't leave the exon.
    if (end > start) {
        int idx = overlapping_region_num(start, end, exons);
        worst = most_severe(worst, splice_boundary_cq(start, end, idx));
    }
    return worst;
}

// realign alt against ref sequence, and strip down to alternate bases
void Tx::trim_alleles(int& start, int&end, std::string& ref, std::string& alt) {
    // trim the front of the alleles, leaving at least one ref base
    size_t k = 0;
    while ((k < alt.size()) && (ref.size() - k > 1) && (alt[k] == ref[k])) {
        k++;
    }
    ref.erase(0, k);
    alt.erase(0, k);
    start += (int) k;
    
    // trim the back of the alleles
    while (!alt.empty() && (ref.size() > 1) && (alt.back() == ref.back())) {
        ref.pop_back();
        alt.pop_back();
        end -= 1;
    }
}

// kept out of line, so building the error message doesn't slow consequence()
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, cold))
#endif
void Tx::ref_mismatch(int pos, const std::string & ref, const std::string & genomic) {
    throw std::invalid_argument("ref allele " + ref + " at " + std::to_string(pos)
        + " doesn't match the genomic sequence for " + get_name() + ": " + genomic);
}

std::string Tx::consequence(int pos, std::string ref, std::string alt) {
    int start = pos;
    int end = pos + (int)ref.size() - 1;
    if (ref.empty()) {
        throw std::invalid_argument("reference allele is empty");
    }
    to_upper(ref);
    to_upper(alt);
    
    // check the ref allele matches the genomic sequence, where they overlap
    if (!genomic_sequence.empty()) {
        int lo = std::max(start, tx_start - gdna_offset);
        int hi = std::min(end, tx_end + gdna_offset);
        const char * genomic = genomic_sequence.data() + genomic_index(lo);
        const char * bases = ref.data() + (lo - start);
        for (int i = 0; i <= hi - lo; i++) {
            if (genomic[i] != bases[i]) {
                ref_mismatch(pos, ref, std::string(genomic, hi - lo + 1));
            }
        }
    }
    
    trim_alleles(start, end, ref, alt);
    
    // insertions lie between the anchor base and its neighbour
    if (ref.size() == 1 && alt.size() > 1 && (alt[0] == ref[0] || alt.back() == ref[0])) {
        if (alt[0] == ref[0]) {
            return insertion_cq(start, alt.substr(1));
        }
        return insertion_cq(start - 1, alt.substr(0, alt.size() - 1));
    }
    
    if ((end < tx_start) || (start > tx_end)) {
        return outside_gene_cq(start, end);
    }
    
    if (exons.size() == 0) {
        throw std::invalid_argument("no exons assigned to this transcript");
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
