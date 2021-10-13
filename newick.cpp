//
// Created by 박석환 on 2021/08/09.
//

#include "newick.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>


void update_nodes(newick_graph *nodes, size_t &count, std::string &lines, size_t &left_pos, size_t &right_pos, newick_start *start_point, newick_graph *end_point) {
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
        start_point->connect_start(start_point, nodes, count, end_point);
    }
    if (check_second == count) {
        start_point->connect_start(start_point, nodes, count + 1, end_point);
    }
    count += 2;
    lines.replace(left_pos, right_pos - left_pos + 1, std::to_string(count));
    right_pos = lines.find(')');
}

void process_newick(const char * newick_file_path, newick_start *start_point, newick_graph *end_point, int &species_num, int *order_max) {
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
    int species_tot_num = std::count(lines.begin(), lines.end(), ':');
    newick_graph *nodes = new newick_graph[species_tot_num];
    size_t left_pos;
    size_t right_pos;
    right_pos = lines.find(')');
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
    species_num = start_point->next.size();
    //give the order (or height) of each node
    end_point->order = -1;
    std::vector<newick_graph*> iterate_previous = end_point->previous;
    while (!iterate_previous.empty()) {
        size_t size = iterate_previous.size();
        for (size_t num = 0; num < size; num++) {
            iterate_previous[num]->order = *order_max;
            if (!iterate_previous[num]->previous.empty()) {
                iterate_previous.insert(iterate_previous.end(), iterate_previous[num]->previous.begin(), iterate_previous[num]->previous.end());
            }
        }
        iterate_previous.erase(iterate_previous.begin(), iterate_previous.begin() + size);
        (*order_max)++;
    }
    (*order_max)--;
}