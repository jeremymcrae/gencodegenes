
#include <algorithm>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "gtf.h"

namespace gencode {

const std::string gene_id_key = "gene_id";
const std::string gene_name_key = "gene_name";
const std::string hgnc_id_key = "hgnc_id";

static inline bool is_trim_char(char c) {
    return c == ';' || c == '=' || c == ' ' || c == '"';
}

// trims separator, space and quote characters from either end of a range.
//
// This operates on a selected range of a GTF line, in a portion where we have 
// identified a value to extract. This avoids creating extra string objects and 
// redundantly running substr.
//
// @param s string for a full GTF line (without line-ending though)
// @param start position where the substring starts
// @param end position where the substring ends (exclusive, may be npos)
// @param out string to assign the trimmed range to, reusing its capacity
static void trim(const std::string &s, size_t start, size_t end, std::string &out) {
    end = std::min(end, s.size());
    while (start < end && is_trim_char(s[start])) {
        start++;
    }
    while (end > start && is_trim_char(s[end - 1])) {
        end--;
    }
    out.assign(s, start, end - start);
}

static std::string trim(const std::string &s, size_t start, size_t end) {
    std::string out;
    trim(s, start, end, out);
    return out;
}

// parse the full attributes field into a key/value map
//
// The attributes field is the 9th (final) tab-delimited column of a GTF line. It
// is a list of "key value" pairs, conventionally separated by "; " (though the
// trailing space is not guaranteed). Values may be quoted. Some keys (e.g. "tag")
// can be repeated; repeated values are joined with commas into a single value.
//
// @param line string for a full GTF line (without line-ending)
// @param offset position where the attributes field starts
static std::map<std::string, std::string> parse_attributes(const std::string &line, size_t offset) {
    std::map<std::string, std::string> attributes;

    // the caller's offset can point at the tab preceding the attributes field,
    // so advance past any leading tabs/spaces to the first attribute character
    offset = line.find_first_not_of("\t ", offset);
    if (offset == std::string::npos) {
        return attributes;
    }

    // the attributes field is the final column, running to the end of the line
    // (a trailing tab, if any, terminates it)
    size_t line_end = line.find('\t', offset);
    if (line_end == std::string::npos) {
        line_end = line.size();
    }
    // drop any trailing line-ending characters from the field
    while (line_end > offset &&
            (line[line_end - 1] == '\n' || line[line_end - 1] == '\r')) {
        line_end--;
    }

    size_t pos = offset;
    while (pos < line_end) {
        // each pair is terminated by a ';' (or the end of the field)
        size_t sep = line.find(';', pos);
        if (sep == std::string::npos || sep > line_end) {
            sep = line_end;
        }

        // the key is the first whitespace-delimited token in [pos, sep)
        size_t key_start = line.find_first_not_of(" ", pos);
        if (key_start == std::string::npos || key_start >= sep) {
            pos = sep + 1;
            continue;
        }
        size_t key_end = line.find(' ', key_start);
        if (key_end == std::string::npos || key_end > sep) {
            // a bare key with no value: keep it as a key with an empty value, so
            // downstream code can probe for the presence of such keys
            std::string bare_key = trim(line, key_start, sep);
            if (!bare_key.empty()) {
                attributes.emplace(std::move(bare_key), "");
            }
            pos = sep + 1;
            continue;
        }

        std::string key = line.substr(key_start, key_end - key_start);
        std::string value = trim(line, key_end, sep);

        auto it = attributes.find(key);
        if (it == attributes.end()) {
            attributes.emplace(std::move(key), std::move(value));
        } else if (!value.empty()) {
            // handle repeated keys (e.g. tag) by joining values with commas
            it->second += "," + value;
        }

        pos = sep + 1;
    }

    return attributes;
}

// check if a whole attribute key matches at a position, returning the position
// after the key, or npos if it doesn't match
static size_t match_key(const std::string &line, size_t pos, const char *key, size_t len) {
    size_t end = pos + len;
    if (end < line.size() && std::memcmp(line.data() + pos, key, len) == 0
            && (line[end] == ' ' || line[end] == '=')) {
        return end;
    }
    return std::string::npos;
}

// find the transcript_id and transcript_type values in one pass over the
// attributes field, matching whole keys only (so "transcript_id" doesn't match
// "havana_transcript_id", or "transcript_type" inside a value). The type falls
// back to transcript_biotype (as in non-gencode GTFs) if transcript_type is
// absent or empty. Fields are set to empty strings if their keys are absent.
static void find_tx_fields(const std::string &line, size_t offset, GTFLine &info) {
    const size_t npos = std::string::npos;
    const char *data = line.data();
    const size_t size = line.size();
    const char prefix[] = "transcript_";
    const size_t prefix_len = sizeof(prefix) - 1;
    bool tx_found = false, type_found = false;
    size_t biotype_end = npos;

    size_t pos = offset;
    while (pos < size) {
        const char *p = static_cast<const char *>(std::memchr(data + pos, 't', size - pos));
        if (p == nullptr) {
            break;
        }
        pos = p - data;
        bool at_start = pos == offset || data[pos - 1] == ' ' || data[pos - 1] == ';' || data[pos - 1] == '\t';
        if (!at_start || pos + prefix_len >= size || data[pos + 1] != 'r'
                || std::memcmp(p, prefix, prefix_len) != 0) {
            pos += 1;
            continue;
        }
        size_t key_pos = pos + prefix_len;
        size_t end;
        if (!tx_found && (end = match_key(line, key_pos, "id", 2)) != npos) {
            trim(line, end, line.find(';', end), info.tx_id);
            tx_found = true;
        } else if (!type_found && (end = match_key(line, key_pos, "type", 4)) != npos) {
            trim(line, end, line.find(';', end), info.transcript_type);
            type_found = true;
        } else if (biotype_end == npos && (end = match_key(line, key_pos, "biotype", 7)) != npos) {
            biotype_end = end;
        } else {
            end = key_pos;
        }
        if (tx_found && type_found && !info.transcript_type.empty()) {
            return;
        }
        pos = end;
    }

    if (!tx_found) {
        info.tx_id.clear();
    }
    if (!type_found) {
        info.transcript_type.clear();
    }
    if (info.transcript_type.empty() && biotype_end != npos) {
        trim(line, biotype_end, line.find(';', biotype_end), info.transcript_type);
    }
}

// parse the required fields from the attributes field
//
// @param all_fields whether to parse the gene fields and attributes map on lines
//     other than "transcript" lines
static void get_attributes_fields(GTFLine &info, std::string &line, int offset, bool all_fields) {
    // tx_id and transcript_type are read for every permitted GTF line in
    // load_transcripts (for transcript-boundary detection and the coding-only
    // filter respectively), so they must be extracted on every line.
    find_tx_fields(line, offset, info);

    // clear fields left over from a previous line, as GTFLine objects are reused
    info.symbol.clear();
    info.alternate_ids.clear();
    info.is_canonical = 0;
    info.attributes.clear();

    // The remaining fields (gene symbol, gene_id, hgnc_id and canonical status)
    // are only consumed by load_transcripts from the "transcript" feature line
    // (or the first line of a transcript, for GTFs without transcript lines),
    // so we defer that work to those lines to avoid redundant string scanning on
    // the many other feature lines. They are derived from the parsed attribute
    // map, rather than re-scanning the line, to avoid extracting fields twice.
    if (all_fields || info.feature == "transcript") {
        info.attributes = parse_attributes(line, offset);

        auto gene_name_it = info.attributes.find(gene_name_key);
        if (gene_name_it != info.attributes.end()) {
            info.symbol = gene_name_it->second;
        }

        auto gene_id_it = info.attributes.find(gene_id_key);
        if (gene_id_it != info.attributes.end() && gene_id_it->second.size() > 0) {
            if (info.symbol.size() == 0) {
                // if we don't have a gene symbol available, use the gene_id field.
                // This means the gene.symbol in the python code can hold non-HGNC
                // data, but it's better than the transcript having a blank name and
                // being used in a gene object with all other blank name transcripts
                info.symbol = gene_id_it->second;
            } else {
                info.alternate_ids.push_back(gene_id_it->second);
            }
        }

        auto hgnc_id_it = info.attributes.find(hgnc_id_key);
        if (hgnc_id_it != info.attributes.end() && hgnc_id_it->second.size() > 0) {
            info.alternate_ids.push_back(hgnc_id_it->second);
        }

        // canonical status is flagged via the "tag" attribute (repeated tag
        // values are joined into a single comma-separated string by parse_attributes)
        auto tag_it = info.attributes.find("tag");
        if (tag_it != info.attributes.end()) {
            if (tag_it->second.find("Ensembl_canonical") != std::string::npos) {
                info.is_canonical = 10;
            } else if (tag_it->second.find("appris_principal") != std::string::npos) {
                info.is_canonical = 5;
            }
        }
    }
}

// parse an integer field from a GTF line, without allocating a substring
//
// @param line GTF line
// @param start position where the field starts
// @param end position where the field ends (the following tab)
static int parse_int(const std::string &line, size_t start, size_t end) {
    // fast path for plain digits, which strtol is slow to handle
    if (end > start && end - start <= 9) {
        int value = 0;
        size_t i = start;
        for (; i < end && line[i] >= '0' && line[i] <= '9'; i++) {
            value = value * 10 + (line[i] - '0');
        }
        if (i == end) {
            return value;
        }
    }

    const char * first = line.c_str() + start;
    char * last;
    errno = 0;
    long value = std::strtol(first, &last, 10);
    if (last == first || last > line.c_str() + end) {
        throw std::invalid_argument("invalid position in GTF line: " + line);
    }
    if (errno == ERANGE || value > INT_MAX || value < INT_MIN) {
        throw std::out_of_range("position out of range in GTF line: " + line);
    }
    return (int) value;
}

// parse required fields from a GTF line
//
// @param line GTF line (without line ending)
// @param all_fields whether to parse the gene fields and attributes map, even
//     if the line isn't a "transcript" line
// @param info GTFLine to fill, reusing the capacity of its strings
void parse_gtfline(std::string & line, GTFLine & info, bool all_fields) {
    // find the tabs ending each of the first 8 fields (chrom, source, feature,
    // start, end, score, strand, frame). The attributes field runs from the
    // final tab to the end of the line. Searching for tabs and extracting
    // substrings is quicker than getline() with a tab delimiter (by 2X).
    size_t tabs[8];
    size_t pos = 0;
    for (int i = 0; i < 8; i++) {
        pos = line.find('\t', pos);
        if (pos == std::string::npos) {
            throw std::invalid_argument("GTF line has fewer than 9 fields: " + line);
        }
        tabs[i] = pos;
        pos += 1;
    }

    info.chrom.assign(line, 0, tabs[0]);
    info.feature.assign(line, tabs[1] + 1, tabs[2] - tabs[1] - 1);
    info.start = parse_int(line, tabs[2] + 1, tabs[3]);
    info.end = parse_int(line, tabs[3] + 1, tabs[4]);
    info.strand.assign(line, tabs[5] + 1, tabs[6] - tabs[5] - 1);

    get_attributes_fields(info, line, tabs[7] + 1, all_fields);
}

GTFLine parse_gtfline(std::string & line, bool all_fields) {
    GTFLine info;
    parse_gtfline(line, info, all_fields);
    return info;
}

GzReader::GzReader(const std::string &path) : file(zng_gzopen(path.c_str(), "rb")), path(path) {}

GzReader::~GzReader() {
    if (file) {
        zng_gzclose(file);
    }
}

// read the next line, without the trailing newline
//
// @param line string to fill with the line contents
// @returns false once the end of the file is reached
bool GzReader::getline(std::string &line) {
    line.clear();
    while (true) {
        if (buf_pos >= buf_len) {
            buf_len = zng_gzread(file, buffer, sizeof(buffer));
            buf_pos = 0;
            if (buf_len <= 0) {
                int32_t errnum = 0;
                const char *msg = zng_gzerror(file, &errnum);
                if (buf_len < 0 || (errnum != Z_OK && errnum != Z_STREAM_END)) {
                    throw std::invalid_argument(msg && *msg ? msg : "error reading GTF: " + path);
                }
                return !line.empty();
            }
        }
        char *newline = static_cast<char *>(std::memchr(buffer + buf_pos, '\n', buf_len - buf_pos));
        if (newline != nullptr) {
            int len = newline - (buffer + buf_pos);
            line.append(buffer + buf_pos, len);
            buf_pos += len + 1;
            return true;
        }
        line.append(buffer + buf_pos, buf_len - buf_pos);
        buf_pos = buf_len;
    }
}

// open GTF file handle
GTF::GTF(std::string path) : reader(path) {
    if (!reader.is_open()) {
        throw std::invalid_argument("cannot open GTF: " + path);
    }
}

// get the next feature line from the GTF, skipping comments and blank lines
//
// @param info GTFLine to fill with the parsed line
// @returns false once the end of the file is reached
bool GTF::next(GTFLine &info) {
    while (reader.getline(line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.find_first_not_of(" \t\r") == std::string::npos || line[0] == '#') {
            continue;
        }
        parse_gtfline(line, info, false);
        if (info.tx_id != prev_tx_id && info.feature != "transcript" && !info.tx_id.empty()) {
            // GTFs without transcript lines need the gene fields from the
            // first line of each transcript
            parse_gtfline(line, info, true);
        }
        prev_tx_id = info.tx_id;
        return true;
    }
    return false;
}

} // namespace
