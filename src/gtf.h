#ifndef GENCODEGENES_GTF_H_
#define GENCODEGENES_GTF_H_

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include <zlib.h>

namespace gencode {

// read lines from a file via zlib, which handles gzipped or uncompressed files
class GzReader {
    gzFile file;
    std::string path;
    char buffer[65536];
    int buf_len = 0;
    int buf_pos = 0;
public:
    explicit GzReader(const std::string &path);
    ~GzReader();
    GzReader(const GzReader &) = delete;
    GzReader & operator=(const GzReader &) = delete;
    bool is_open() const { return file != nullptr; }
    bool getline(std::string &line);
};

// store required fields from a GTF line
struct GTFLine {
    std::string chrom;
    std::string feature;
    int start;
    int end;
    std::string strand;
    std::string symbol;
    std::vector<std::string> alternate_ids;
    std::string tx_id;
    std::string transcript_type;
    int is_canonical = 0;
    std::map<std::string, std::string> attributes;
};

GTFLine parse_gtfline(std::string &line, bool all_fields=false);
void parse_gtfline(std::string &line, GTFLine &info, bool all_fields);

class GTF
{
    GzReader reader;
    std::string line;
    std::string prev_tx_id;
public:
    GTF(std::string path);
    bool next(GTFLine &info);
};

} // namespace

#endif // GENCODEGENES_GTF_H_
