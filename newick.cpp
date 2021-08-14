//
// Created by 박석환 on 2021/08/09.
//

#include "newick.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>


void update_nodes(newick_graph *nodes, size_t &count, std::string &lines, size_t &left_pos, size_t &right_pos, newick_start start_point, newick_graph end_point) {
    size_t temp_pos_comma;
    size_t temp_pos_colon_first;
    size_t temp_pos_colon_second;
    std::string first_species;
    std::string second_species;
    std::string first_bl_string;
    std::string second_bl_string;
    float first_bl;
    float second_bl;
    std::string parent;
    temp_pos_comma = lines.find(',', left_pos);
    temp_pos_colon_first = lines.find(':', left_pos);
    temp_pos_colon_second = lines.find(':', temp_pos_comma);
    first_species = lines.substr(left_pos + 1, temp_pos_colon_first - left_pos - 1);
    second_species = lines.substr(temp_pos_comma + 1, temp_pos_colon_second - temp_pos_comma - 1);
    first_bl_string = lines.substr(temp_pos_colon_first + 1, temp_pos_comma - temp_pos_colon_first - 1);
    first_bl = std::stof(first_bl_string);
    second_bl_string = lines.substr(temp_pos_colon_second + 1, right_pos - temp_pos_colon_second - 1);
    second_bl = std::stof(second_bl_string);
    parent = std::to_string(count + 2);
    nodes[count].parent = parent;
    nodes[count].branch_length = first_bl;
    nodes[count].species = first_species;
    nodes[count + 1].parent = parent;
    nodes[count + 1].branch_length = second_bl;
    nodes[count + 1].species = second_species;
    size_t check_first = 0;
    size_t check_second = 0;
    for (size_t i = 0; i < count; i++) {
        if (nodes[i].parent == nodes[count].species) {
            nodes->insert_inbetween_end_n_target(end_point, nodes, count, nodes[i]);
        } else {
            check_first++;
        }
        if (nodes[i].parent == nodes[count + 1].species) {
            nodes->insert_inbetween_end_n_target(end_point, nodes, count + 1, nodes[i]);
        } else {
            check_second++;
        }
    }
    if (check_first == count) {
        start_point.connect_start(start_point, nodes, count, end_point);
    }
    if (check_second == count) {
        start_point.connect_start(start_point, nodes, count + 1, end_point);
    }
    count += 2;
    lines.replace(left_pos, right_pos - left_pos + 1, std::to_string(count));
    right_pos = lines.find(')');
}

newick_graph *process_newick(const char * newick_file_path) {
    newick_start start_point;
    newick_graph end_point;
    std::ifstream newick;
    newick.open(newick_file_path);
    std::string line;
    std::string lines;
    getline(newick, line);
    while (!line.empty()) {
        lines += line;
        // clean lines (get rid of all the whitespaces)
        while (lines.find(' ') != std::string::npos) {
            size_t pos = lines.find(' ');
            lines.erase(lines.begin() + pos);
        }
        getline(newick, line);
    }
    // make binary tree
    uint16_t nodes_count = std::count(lines.begin(), lines.end(), ':') + 1;
    newick_graph *nodes = new newick_graph[nodes_count];
    size_t left_pos;
    size_t right_pos;
    const char left_delimiter = '(';
    const char right_delimiter = ')';
    right_pos = lines.find(right_delimiter);
    size_t count = 0;
    while (right_pos != std::string::npos) {
        for (size_t pos = right_pos; pos >= 0; pos--) {
            if (lines[pos] == '(') {
                left_pos = pos;
                break;
            }
        }
        update_nodes(nodes, count, lines, left_pos, right_pos, start_point, end_point);
    }


    return nodes;
}