//
// Created by 박석환 on 2021/08/15.
//

#include "process_maf.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
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
                intact_codon = true;
                if (positive[start + num] == false) {
                    continue;
                } else if (positive[start + num] == true) {
                    if (positive[start + num + 1] == true && positive[start + num + 2] == true) {
                        codon = block[0].substr(num, 3);
                        std::transform(codon.begin(), codon.end(), codon.begin(), ::toupper);
                        if (codon_to_char.find(codon) == codon_to_char.end()) {
                            num += 2;
                            continue;
                        } else {
                            codon_set[0].codon = codon_to_char.find(codon)->second;
                            for (size_t block_it = 1; block_it < species_num; block_it++) {
                                codon = block[block_it].substr(num, 3);
                                std::transform(codon.begin(), codon.end(), codon.begin(), ::toupper);
                                if (codon_to_char.find(codon) == codon_to_char.end()) {
                                    intact_codon = false;
                                    break;
                                } else {
                                    codon_set[block_it].codon = codon_to_char.find(codon)->second;
                                }
                            }
                        }
                        if (intact_codon == true) {
                            coding_codon_set.emplace_back(codon_set);
                        }
                        num += 2;
                    }
                } else {
                    throw ("not being bool in vector of bool");
                }
            }
            for (size_t num = 0; num < len - 2; num++) {
                intact_codon = true;
                if (negative[start + num] == false) {
                    continue;
                } else if (negative[start + num] == true) {
                    if (negative[start + num + 1] == true && negative[start + num + 2] == true) {
                        codon = block[0].substr(num, 3);
                        std::transform(codon.begin(), codon.end(), codon.begin(), ::toupper);
                        return_complement_codon(codon);
                        if (codon_to_char.find(codon) == codon_to_char.end()) {
                            num += 2;
                            continue;
                        } else {
                            codon_set[0].codon = codon_to_char.find(codon)->second;
                            for (size_t block_it = 1; block_it < species_num; block_it++) {
                                codon = block[block_it].substr(num, 3);
                                std::transform(codon.begin(), codon.end(), codon.begin(), ::toupper);
                                return_complement_codon(codon);
                                if (codon_to_char.find(codon) == codon_to_char.end()) {
                                    intact_codon = false;
                                    break;
                                } else {
                                    codon_set[block_it].codon = codon_to_char.find(codon)->second;
                                }
                            }
                        }
                        if (intact_codon == true) {
                            coding_codon_set.emplace_back(codon_set);
                        }
                        num += 2;
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