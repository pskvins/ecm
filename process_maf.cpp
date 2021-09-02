//
// Created by 박석환 on 2021/08/15.
//

#include "process_maf.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include "codon_frequency.cpp"

std::vector<aligned_codon> maf_to_aligned_codons(const char * path_to_maf, short int species_num) {
    std::ifstream maf;
    maf.open(path_to_maf);
    std::string line;
    std::vector<std::string> block;
    aligned_codon *codon_set = new aligned_codon[species_num];
    std::vector<aligned_codon> codon_set_vector;
    size_t pos = 0;
    size_t pos_2 = 0;
    std::string codon;
    char codon_map;
    getline(maf, line);
    while (!line.empty()) {
        if (line[0] == 'a') {
            getline(maf, line);
            while (line[0] != 'a') {
                if (line [0] == 's') {
                    block.emplace_back(line);
                }
            }
            for (size_t num = 0; num < block.size(); num++) {
                pos = block[0].find('\t');
                pos_2 = block[0].find('.');
                codon_set[0].species = block[0].substr(pos + 1, pos_2 - pos - 1);
                pos = block[0].find_last_of('\t');
                block[0].erase(0, pos);
            }
            while (block[0].size() >= 3) {
                for (size_t num = 0; num < block.size(); num++) {
                    codon = block[num].substr(0, 3);
                    codon_map = codon_to_char.find(codon)->second;
                    codon_set[num].codon = codon_map;
                    block[num].erase(0, 3);
                }
                codon_set_vector.emplace_back(*codon_set);
            }
        } else {
            while (line[0] != 'a') {
                getline(maf, line);
            }
        }
    }
    return codon_set_vector;
}