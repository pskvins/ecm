#include "njtree.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>


uint32_t finding_index(std::string * name,std::string chrom) {
    int pos_chrom = 0;
    for (size_t pos = 0; pos < sizeof(name); pos++) {
        if (name[pos] == ">" + chrom) {
            break;
        }
        else {
            pos_chrom++;
        }
    }
    return pos_chrom;
}

struct maf_info {
    std::string name;
    std::string seq;
    uint32_t length;

    maf_info (std::string name, std::string seq, uint32_t length) :
        name(name),
        seq(seq),
        length(length) {}
};

void read_maf(std::string delimiter, std::string line, std::vector<maf_info> &block) {
    size_t pos;
    uint queue;
    std::string name;
    std::string seq;
    uint32_t length;
    while (queue < 7) {
        pos = line.find(delimiter);
        if (queue == 1) {
            name = line.substr(0, line.find("."));
        } else if (queue == 6) {
            seq = line.substr(0, pos);
            length = seq.size();
        }
        queue++;
        line.erase(0, pos + delimiter.length());
    }
    block.emplace_back(name, seq, length);
}

uint32_t ** from_maf_to_distance_matrix (const char *file_path, int species_num) {
    uint32_t **distance_matrix = new uint32_t * [species_num];
    for (size_t i = 0; i < species_num; i++) {
        distance_matrix[i] = new uint32_t [species_num];
    }
    std::string names[species_num];
    for (size_t i = 0; i < sizeof(names); i++) {
        names[i] = "Null";
    }
    std::ifstream maf;
    maf.open(file_path);
    std::string line;
    std::vector<maf_info> block;
    std::vector<std::string> spe_names;
    std::string delimiter = "\t";
    getline(maf, line);
    if (line[0] == '#') {
        while (line[0] == '#') {
            getline(maf, line);
        }
    }
    while (!(line.empty())) {
        if (line[0] == '#') {
            while (line[0] == '#') {
                getline(maf, line);
            }
        }
        if (line[0] == 'a') {
            block.clear();
            getline(maf, line);
            while (line[0] == 's') {
                read_maf(delimiter, line, block);
                getline(maf, line);
            }
            if (std::find(names, names + species_num, "Null") != names + species_num) {
                for (size_t i = 0; i < block.size() ; i++) {
                    if (std::find(names, names + species_num, block[i].name) == names + species_num) {
                        for (size_t j = 0; j < sizeof(names); j++) {
                            if (names[j] == "Null") {
                                names[j] = block[i].name;
                                break;
                            }
                        }
                    }
                }
            }
            for (size_t row = 0; row < block.size(); row++) {
                for (size_t col = 0; col < block.size(); col++) {
                    if (row < col) {
                        assert(block[row].length == block[col].length);
                        for (size_t pos = 0; pos < block[row].length; pos++) {
                            if (block[row].seq[pos] != block[col].seq[pos] && block[row].seq[pos] != 'N' && block[row].seq[pos] != 'n'
                             && block[row].seq[pos] != '-' && block[col].seq[pos] != 'N' && block[col].seq[pos] != 'n' && block[col].seq[pos] != '-') {
                                ++distance_matrix[finding_index(names, block[row].name)][finding_index(names, block[col].name)];
                            }
                        }
                    }
                }
            }
        } else {
            getline(maf, line);
        }
    }
    for (size_t i = 0; i < species_num; i++) {
        delete [] distance_matrix[i];
    }
    delete [] distance_matrix;
    return distance_matrix;
}

