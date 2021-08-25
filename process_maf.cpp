//
// Created by 박석환 on 2021/08/15.
//

#include "process_maf.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include "codon_frequency.cpp"

void process_reference(std::string &line, &all_CDS_info) {
    uint32_t start;
    uint32_t end;
    size_t pos_start;
    size_t pos_end;
    char delimiter = '\t';

}


std::vector<aligned_codon> maf_to_aligned_codons(const char * path_to_maf, std::vector<CDS_info> &all_CDS_info) {
    std::ifstream maf;
    maf.open(path_to_maf);
    std::string line;
    std::vector<std::string> block;
    getline(maf, line);
    while (!line.empty()) {
        if (line[0] == '#') {
            while (line[0] == '#') {
                getline(maf, line);
            }
        }
        if (line[0] == 'a') {
            getline(maf, line);

            while (line[0] == 's') {

            }
        }
    }
}