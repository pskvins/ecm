//
// Created by 박석환 on 2021/08/09.
//

#include "expectation_step.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <iostream>

newick_start start_point();
newick_graph end_point();

void process_newick(const char * newick_file_path) {
    std::ifstream newick;
    newick.open(newick_file_path);
    std::string line;
    std::string lines;
    getline(newick, line);
    while (!line.empty()) {
        lines += line;
        getline(newick, line);
    }
    uint16_t nodes_count = std::count(lines.begin(), lines.end(), ':') + 1;
    newick_graph *nodes = new newick_graph[nodes_count];
    size_t temp_pos;
    size_t left_pos;
    size_t right_pos;
    const char left_delimiter = '(';
    const char right_delimiter = ')';
    const char delimit_name_bl = ':';
    std::vector<size_t> left_delimit_pos;
    left_pos = lines.find(left_delimiter);
    left_delimit_pos.emplace_back(left_pos);
    right_pos = lines.find(right_delimiter);
    size_t count = 0;
    while (left_pos != std::string::npos) {
        left_pos = lines.find(left_delimiter, left_pos + 1);
        left_delimit_pos.emplace_back(left_pos);
        if (left_delimit_pos.back() > right_pos) {
            while (left_delimit_pos.back() > right_pos) {

            }
        }
    }
}

float possibility_change_codon(std::string codon_pre, std::string codon_pro, float branch_length, float **q_matrix) {

}