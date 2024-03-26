
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "gzstream/gzstream.h"

#include "gtf.h"

namespace gencode {

std::string trim(const std::string &s, const std::string &vals) {
    size_t start = s.find_first_not_of(vals);
    if (start == std::string::npos) {
        return "";
    }
    size_t end = s.find_last_not_of(vals);
    return s.substr(start, (end + 1) - start);
}

// parse the required fields from the attributes field
static void get_attributes_fields(GTFLine &info, std::string &line, int offset) {
    const std::string tx_id_key = "transcript_id ";
    const std::string gene_id_key = "gene_id ";
    const std::string gene_name_key = "gene_name ";
    std::string type_key = "transcript_type ";
    const std::string hgnc_id_key = "hgnc_id ";

    size_t tx_start = line.find(tx_id_key, offset) + tx_id_key.size();
    size_t tx_end = line.find(";", tx_start);

    if (tx_start - tx_id_key.size() == std::string::npos) {
        // handle if the string was not found
        tx_start = offset;
        tx_end = offset;
    }

    size_t gene_id_start = line.find(gene_id_key, offset) + gene_id_key.size();
    size_t gene_id_end = line.find(";", gene_id_start);
    if (gene_id_start - gene_id_key.size() == std::string::npos) {
        // handle if the string was not found
        gene_id_start = offset;
        gene_id_end = offset;
    }
    
    size_t gene_start = line.find(gene_name_key, offset) + gene_name_key.size();
    size_t gene_end = line.find(";", gene_start);

    if (gene_start - gene_name_key.size() == std::string::npos) {
        // handle if the string was not found
        gene_start = tx_end;
        gene_end = tx_end;
    }

    size_t type_start = line.find(type_key, offset) + type_key.size();
    if (type_start - type_key.size() == std::string::npos) {
        // allow for alternate transcript_type key, as found in non-gencode GTF files 
        type_key = "transcript_biotype ";
        type_start = line.find(type_key, offset) + type_key.size();
    }
    size_t type_end = line.find(";", type_start);

    if (type_start - type_key.size() == std::string::npos) {
        // handle if the string was not found
        type_start = gene_end;
        type_end = gene_end;
    }

    size_t hgnc_id_start = line.find(hgnc_id_key, offset) + hgnc_id_key.size();
    size_t hgnc_id_end = line.find(";", hgnc_id_start);
    if (hgnc_id_start - hgnc_id_key.size() == std::string::npos) {
        // handle if the string was not found
        hgnc_id_start = offset;
        hgnc_id_end = offset;
    }
    
    info.symbol = trim(line.substr(gene_start, gene_end - gene_start), " \"");
    info.tx_id = trim(line.substr(tx_start, tx_end - tx_start), " \"");
    info.transcript_type = trim(line.substr(type_start, type_end - type_start), " \"");
    
    if (gene_id_start != gene_id_end) {
        std::string gene_id = trim(line.substr(gene_id_start, gene_id_end - gene_id_start), " \"");
        if (info.symbol.size() == 0) {
            info.symbol = gene_id;
        } else {
            info.alternate_ids.push_back(gene_id);
        }
    }
    if (hgnc_id_start != hgnc_id_end) {
        info.alternate_ids.push_back(trim(line.substr(hgnc_id_start, hgnc_id_end - hgnc_id_start), " \""));
    }
    
    if (info.feature == "transcript") {
        if (line.find("appris_principal", offset) != std::string::npos) {
            info.is_canonical = 5;
        } else if (line.find("Ensembl_canonical", offset) != std::string::npos) {
            info.is_canonical = 10;
        }
    }
}

// parse required fields from a GTF line
GTFLine parse_gtfline(std::string & line) {
    if (line.size() == 0) {
        throw std::out_of_range("end of file");
    }

    GTFLine info;

    // there are only a few fields we need from the GTF lines, and some fields
    // are only a single character long, so it's quickest to search for the next
    // tab along the line, then extract the substring to get the required fields.
    // getline() with tab delimiter was 2X slower.
    int chr_idx = 0;
    int source_idx = line.find("\t", chr_idx);
    int feature_idx = line.find("\t", source_idx + 6);
    int start_idx = line.find("\t", feature_idx + 3);
    int end_idx = line.find("\t", start_idx + 2);
    int score_idx = line.find("\t", end_idx + (end_idx - start_idx));

    info.chrom = line.substr(chr_idx, source_idx - chr_idx);
    info.feature = line.substr(feature_idx + 1, start_idx - feature_idx - 1);
    info.start = std::stoi(line.substr(start_idx + 1, end_idx - start_idx - 1));
    info.end = std::stoi(line.substr(end_idx + 1, score_idx - end_idx - 1));
    info.strand = line[score_idx + 3];

    get_attributes_fields(info, line, score_idx + 6);

    return info;
}

// open GTF file handle
GTF::GTF(std::string path) {
    gzipped = path.substr(path.length()-2, 2) == "gz";
    if (gzipped) {
        gzhandle.open(path.c_str());
    } else {
        handle.open(path, std::ios::in);
    }
}

// get next line from the GTF
GTFLine GTF::next() {
    if (gzipped) { 
        std::getline(gzhandle, line);
    } else { 
        std::getline(handle, line);
    }
    while (line[0] == '#') {
        if (gzipped) {
            std::getline(gzhandle, line);
        } else {
            std::getline(handle, line);
        }
    }
    return parse_gtfline(line);
}

} // namespace
