
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

const std::string tx_id_key = "transcript_id";
const std::string type_key = "transcript_type";
const std::string biotype_key = "transcript_biotype";
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
static std::string trim(const std::string &s, size_t start, size_t end) {
    end = std::min(end, s.size());
    while (start < end && is_trim_char(s[start])) {
        start++;
    }
    while (end > start && is_trim_char(s[end - 1])) {
        end--;
    }
    return s.substr(start, end - start);
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

// find the value for an attribute key, matching whole keys only (so "transcript_id"
// doesn't match "havana_transcript_id", or "transcript_type" inside a value)
//
// @returns the value, or an empty string if the key is absent
static std::string find_attribute(const std::string &line, const std::string &key, size_t offset) {
    size_t pos = offset;
    while ((pos = line.find(key, pos)) != std::string::npos) {
        size_t end = pos + key.size();
        bool at_start = pos == offset || line[pos - 1] == ' ' || line[pos - 1] == ';' || line[pos - 1] == '\t';
        if (at_start && end < line.size() && (line[end] == ' ' || line[end] == '=')) {
            return trim(line, end, line.find(';', end));
        }
        pos = end;
    }
    return "";
}

// parse the required fields from the attributes field
//
// @param all_fields whether to parse the gene fields and attributes map on lines
//     other than "transcript" lines
static void get_attributes_fields(GTFLine &info, std::string &line, int offset, bool all_fields) {
    // tx_id and transcript_type are read for every permitted GTF line in
    // load_transcripts (for transcript-boundary detection and the coding-only
    // filter respectively), so they must be extracted on every line.
    info.tx_id = find_attribute(line, tx_id_key, offset);
    info.transcript_type = find_attribute(line, type_key, offset);
    if (info.transcript_type.empty()) {
        // allow for alternate transcript_type key, as found in non-gencode GTF files
        info.transcript_type = find_attribute(line, biotype_key, offset);
    }

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
GTFLine parse_gtfline(std::string & line, bool all_fields) {
    GTFLine info;

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

    info.chrom = line.substr(0, tabs[0]);
    info.feature = line.substr(tabs[1] + 1, tabs[2] - tabs[1] - 1);
    info.start = parse_int(line, tabs[2] + 1, tabs[3]);
    info.end = parse_int(line, tabs[3] + 1, tabs[4]);
    info.strand = line.substr(tabs[5] + 1, tabs[6] - tabs[5] - 1);

    get_attributes_fields(info, line, tabs[7] + 1, all_fields);

    return info;
}

GzReader::GzReader(const std::string &path) : file(gzopen(path.c_str(), "rb")), path(path) {}

GzReader::~GzReader() {
    if (file) {
        gzclose(file);
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
            buf_len = gzread(file, buffer, sizeof(buffer));
            buf_pos = 0;
            if (buf_len <= 0) {
                int errnum = 0;
                const char *msg = gzerror(file, &errnum);
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
        info = parse_gtfline(line);
        if (info.tx_id != prev_tx_id && info.feature != "transcript" && !info.tx_id.empty()) {
            // GTFs without transcript lines need the gene fields from the
            // first line of each transcript
            info = parse_gtfline(line, true);
        }
        prev_tx_id = info.tx_id;
        return true;
    }
    return false;
}

} // namespace
