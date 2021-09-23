//
// Created by 박석환 on 2021/08/15.
//

#include "process_maf.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include "codon_frequency.cpp"
#define nuc2int(x) (x & 14u)>>1u

std::vector<std::vector<aligned_codon>> make_msa_to_aligned_coding_codon(const char *file_path_maf, std::vector<CDS_info> coding_regions, int species_num) {

    std::vector<bool> positive;
    positive.reserve(25000000);
    std::vector<bool> negative;
    negative.reserve(25000000);
    std::vector<std::vector<aligned_codon>> coding_codon_set;
    std::vector<CDS_info> coding_regions_chrom;
    std::vector<aligned_codon> codon_set;
    aligned_codon temp;
    for (int i = 0; i < species_num; i++) {
        codon_set.emplace_back(temp);
    }
    std::ifstream maf;
    maf.open(file_path_maf);
    std::string line;
    std::vector<std::string> block;
    int start;
    int len;
    std::string chrom;
    std::string chrom_old = "nothing";
    size_t pos_start;
    size_t pos_middle;
    size_t pos_end;
    size_t pos_dot;
    getline(maf, line);
    std::string codon;
    bool intact_codon;
    while (!line.empty()) {
        if (line[0] == 'a') {
            block.clear();
            getline(maf, line);
            pos_start = line.find('\t');
            pos_dot = line.find('.');
            codon_set[0].species = line.substr(pos_start + 1, pos_dot - pos_start - 1);
            pos_start = line.find('\t', pos_start + 1);
            pos_middle = line.find('\t', pos_start + 1);
            pos_end = line.find('\t', pos_middle + 1);
            pos_dot = line.find_last_of('.', pos_start);
            start = std::stoi(line.substr(pos_start + 1, pos_middle - pos_start - 1));
            len = std::stoi(line.substr(pos_middle + 1, pos_end - pos_middle - 1));
            chrom = line.substr(pos_dot + 1, pos_start - pos_dot - 1);
            pos_start = line.find_last_of('\t');
            line.erase(0, pos_start + 1);
            if (line.size() < 3) {
                continue;
            }
            block.emplace_back(line);
            getline(maf, line);
            if (chrom != chrom_old) {
                coding_regions_chrom.clear();
                std::fill(positive.begin(), positive.end(), false);
                std::fill(negative.begin(), negative.end(), false);
                for (size_t num = 0; num < coding_regions.size(); num++) {
                    if (coding_regions[num].chrom == chrom) {
                        coding_regions_chrom.emplace_back(coding_regions[num]);
                    }
                }
                for (size_t num = 0; num < coding_regions_chrom.size(); num++) {
                    if (coding_regions_chrom[num].chrom == chrom) {
                        if (coding_regions_chrom[num].strand == "+") {
                            for (int start = coding_regions_chrom[num].start + coding_regions_chrom[num].additional; start < coding_regions_chrom[num].end + 1; start++) {
                                positive[start] = true;
                            }
                        } else if (coding_regions_chrom[num].strand == "-") {
                            for (int start = coding_regions_chrom[num].start + coding_regions_chrom[num].additional; start < coding_regions_chrom[num].end + 1; start++) {
                                negative[start] = true;
                            }
                        } else {
                            throw ("strand not being positive or negative");
                        }
                    }
                }
                chrom_old = chrom;
            }
            size_t it = 1;
            while (line[0] == 's') {
                pos_start = line.find('\t');
                pos_dot = line.find('.');
                codon_set[it].species = line.substr(pos_start + 1, pos_dot - pos_start - 1);
                pos_start = line.find_last_of('\t');
                line.erase(0, pos_start + 1);
                block.emplace_back(line);
                getline(maf, line);
                it++;
            }
            for (size_t num = 0; num < len - 2; num++) {
                if (positive[start + num] == false) {
                    continue;
                } else if (positive[start + num] == true) {
                    if (positive[start + num + 1] == true && positive[start + num + 2] == true) {
                        intact_codon = true;
                        codon = block[0].substr(num, 3);
                        std::transform(codon.begin(), codon.end(), codon.begin(), ::toupper);
                        int coordinate[3] = {0, 0, 0};
                        for (short int base = 0; base < 3; base++) {
                            if (codon[base] != 'T' && codon[base] != 'A' && codon[base] != 'C' && codon[base] != 'G'){
                                num += 2;
                                intact_codon = false;
                                break;
                            } else {
                                coordinate[base] = nuc2int(codon[base]);
                            }
                        }
                        if (intact_codon == true) {
                            codon_set[0].codon = codonTable[coordinate[0]][coordinate[1]][coordinate[2]];
                            for (size_t block_it = 1; block_it < species_num && intact_codon == true; block_it++) {
                                codon = block[block_it].substr(num, 3);
                                std::transform(codon.begin(), codon.end(), codon.begin(), ::toupper);
                                for (short int base = 0; base < 3; base++) {
                                    if (codon[base] != 'T' && codon[base] != 'A' && codon[base] != 'C' && codon[base] != 'G') {
                                        num += 2;
                                        intact_codon = false;
                                        break;
                                    } else {
                                        coordinate[base] = nuc2int(codon[base]);
                                    }
                                }
                                codon_set[block_it].codon = codonTable[coordinate[0]][coordinate[1]][coordinate[2]];
                            }
                            if (intact_codon == true) {
                                coding_codon_set.emplace_back(codon_set);
                                num += 2;
                            }
                        }
                    }
                } else {
                    throw ("not being bool in vector of bool");
                }
            }
            for (size_t num = 0; num < len - 2; num++) {
                if (negative[start + num] == false) {
                    continue;
                } else if (negative[start + num] == true) {
                    if (negative[start + num + 1] == true && negative[start + num + 2] == true) {
                        intact_codon = true;
                        codon = block[0].substr(num, 3);
                        std::transform(codon.begin(), codon.end(), codon.begin(), ::toupper);
                        return_complement_codon(codon);
                        int coordinate[3] = {0, 0, 0};
                        for (short int base = 0; base < 3; base++) {
                            if (codon[base] != 'T' && codon[base] != 'A' && codon[base] != 'C' || codon[base] != 'G'){
                                num += 2;
                                intact_codon = false;
                                break;
                            } else {
                                coordinate[base] = nuc2int(codon[base]);
                            }
                        }
                        if (intact_codon == true) {
                            codon_set[0].codon = codonTable[coordinate[0]][coordinate[1]][coordinate[2]];
                            for (size_t block_it = 1; block_it < species_num && intact_codon == true; block_it++) {
                                codon = block[block_it].substr(num, 3);
                                std::transform(codon.begin(), codon.end(), codon.begin(), ::toupper);
                                return_complement_codon(codon);
                                for (short int base = 0; base < 3; base++) {
                                    if (codon[base] != 'T' && codon[base] != 'A' && codon[base] != 'C' && codon[base] != 'G') {
                                        num += 2;
                                        intact_codon = false;
                                        break;
                                    } else {
                                        coordinate[base] = nuc2int(codon[base]);
                                    }
                                }
                                codon_set[block_it].codon = codonTable[coordinate[0]][coordinate[1]][coordinate[2]];
                            }
                            if (intact_codon == true) {
                                coding_codon_set.emplace_back(codon_set);
                                num += 2;
                            }
                        }
                    }
                }
            }
        } else {
            while (line[0] != 'a') {
                getline(maf, line);
            }
        }
    }
    return coding_codon_set;
}