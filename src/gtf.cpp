
#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <map>
#include <stdexcept>
#include <string>

#include "gtf.h"

namespace gencode {

const std::string gene_id_key = "gene_id";
const std::string gene_name_key = "gene_name";
const std::string hgnc_id_key = "hgnc_id";
const std::string tag_key = "tag";

static inline bool is_trim_char(char c) {
    return c == ';' || c == '=' || c == ' ' || c == '"';
}

// assign a range of a GTF line to a string, trimming separator, space and quote
// characters from either end. Assigning to an existing string reuses its capacity.
//
// @param s string for a full GTF line
// @param start position where the range starts
// @param end position where the range ends (exclusive, may be npos)
// @param out string to assign the trimmed range to
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

// parse the full attributes field into a key/value map
//
// The attributes field is the 9th (final) tab-delimited column of a GTF line. It
// is a list of "key value" pairs, conventionally separated by "; " (though the
// trailing space is not guaranteed). Values may be quoted. Some keys (e.g. "tag")
// can be repeated; repeated values are joined with commas into a single value.
// Bare keys without a value are kept with an empty value, so downstream code can
// probe for the presence of such keys.
//
// @param line string for a full GTF line
// @param offset position where the attributes field starts
static std::map<std::string, std::string> parse_attributes(const std::string &line, size_t offset) {
    std::map<std::string, std::string> attributes;

    // the field runs to a trailing tab (if any) or the end of the line, without
    // any line-ending characters
    offset = std::min(line.find_first_not_of("\t ", offset), line.size());
    size_t field_end = std::min(line.find('\t', offset), line.size());
    while (field_end > offset && (line[field_end - 1] == '\n' || line[field_end - 1] == '\r')) {
        field_end--;
    }

    std::string key, value;
    for (size_t pos = offset; pos < field_end; ) {
        size_t sep = std::min(line.find(';', pos), field_end);
        size_t key_start = pos;
        while (key_start < sep && line[key_start] == ' ') {
            key_start++;
        }
        size_t key_end = key_start;
        while (key_end < sep && line[key_end] != ' ') {
            key_end++;
        }
        pos = sep + 1;

        if (key_end < sep) {
            key.assign(line, key_start, key_end - key_start);
            trim(line, key_end, sep, value);
        } else {
            trim(line, key_start, sep, key);
            value.clear();
        }
        if (key.empty()) {
            continue;
        }

        auto it = attributes.find(key);
        if (it == attributes.end()) {
            attributes.emplace(key, value);
        } else if (!value.empty()) {
            it->second += ',';
            it->second += value;
        }
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
//
// These are needed on every line (for transcript boundaries and the coding-only
// filter), so this avoids parsing the full attributes field.
static void find_tx_fields(const std::string &line, size_t offset, GTFLine &info) {
    const size_t npos = std::string::npos;
    const char prefix[] = "transcript_";
    const size_t prefix_len = sizeof(prefix) - 1;
    const char *data = line.data();
    const size_t size = line.size();

    info.tx_id.clear();
    info.transcript_type.clear();
    bool tx_found = false, type_found = false;
    size_t biotype_end = npos;

    for (size_t pos = offset; pos < size; ) {
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
        pos += prefix_len;
        size_t end;
        if (!tx_found && (end = match_key(line, pos, "id", 2)) != npos) {
            trim(line, end, line.find(';', end), info.tx_id);
            tx_found = true;
        } else if (!type_found && (end = match_key(line, pos, "type", 4)) != npos) {
            trim(line, end, line.find(';', end), info.transcript_type);
            type_found = true;
        } else if (biotype_end == npos && (end = match_key(line, pos, "biotype", 7)) != npos) {
            biotype_end = end;
        } else {
            continue;
        }
        if (tx_found && !info.transcript_type.empty()) {
            return;
        }
        pos = end;
    }

    if (info.transcript_type.empty() && biotype_end != npos) {
        trim(line, biotype_end, line.find(';', biotype_end), info.transcript_type);
    }
}

// parse the attributes map, and the gene fields derived from it
//
// These are only used from the first line of each transcript (usually the
// "transcript" line), so are skipped on other lines.
static void parse_gene_fields(const std::string &line, size_t offset, GTFLine &info) {
    info.attributes = parse_attributes(line, offset);
    auto get = [&info](const std::string &key) -> const std::string & {
        static const std::string missing;
        auto it = info.attributes.find(key);
        return it == info.attributes.end() ? missing : it->second;
    };

    // fall back to the gene_id if there is no gene symbol. This means the
    // gene.symbol in the python code can hold non-HGNC data, but it's better than
    // grouping all transcripts lacking a symbol into one blank-named gene.
    info.symbol = get(gene_name_key);
    const std::string &gene_id = get(gene_id_key);
    if (!gene_id.empty()) {
        if (info.symbol.empty()) {
            info.symbol = gene_id;
        } else {
            info.alternate_ids.push_back(gene_id);
        }
    }
    const std::string &hgnc_id = get(hgnc_id_key);
    if (!hgnc_id.empty()) {
        info.alternate_ids.push_back(hgnc_id);
    }

    // canonical status is flagged in the (comma-joined) tag values
    const std::string &tags = get(tag_key);
    if (tags.find("Ensembl_canonical") != std::string::npos) {
        info.is_canonical = 10;
    } else if (tags.find("appris_principal") != std::string::npos) {
        info.is_canonical = 5;
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
// @param info GTFLine to fill, reusing the capacity of its strings
// @param all_fields whether to parse the gene fields and attributes map, even
//     if the line isn't a "transcript" line
void parse_gtfline(std::string & line, GTFLine & info, bool all_fields) {
    // find the tabs ending each of the first 8 fields (chrom, source, feature,
    // start, end, score, strand, frame). The attributes field runs from the
    // final tab to the end of the line.
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
    find_tx_fields(line, tabs[7] + 1, info);

    // clear fields left over from a previous line
    info.symbol.clear();
    info.alternate_ids.clear();
    info.is_canonical = 0;
    info.attributes.clear();
    if (all_fields || info.feature == "transcript") {
        parse_gene_fields(line, tabs[7] + 1, info);
    }
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
