#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>


uint32_t finding_index(std::vector<std::string> &name,std::string chrom) {
    int pos_chrom = 0;
    for (size_t pos = 0; pos < name.size(); pos++) {
        if (name[pos] == chrom) {
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

struct distance_matrix_grid {
    uint32_t row;
    uint32_t col;

    distance_matrix_grid (uint32_t row, uint32_t col) :
        row(row),
        col(col) {}
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

float ** from_maf_to_distance_matrix (std::vector<std::string> &names, const char *file_path, int species_num) {
    float **distance_matrix = new float * [species_num];
    for (size_t i = 0; i < species_num; i++) {
        distance_matrix[i] = new float [species_num];
    }
    for (size_t row = 0; row < species_num; row++) {
        for (size_t col = 0; col < species_num; col++) {
            distance_matrix[row][col] = 0;
        }
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
            if (names.size() != species_num) {
                for (size_t i = 0; i < block.size() ; i++) {
                    if (std::find(names.begin(), names.end(), block[i].name) == names.end()) {
                        names.emplace_back(block[i].name);
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
    // make distance matrix neat
    for (size_t row = 0; row < species_num; row++) {
        for (size_t col = 0; col < species_num; col++) {
            if (row < col) {
                distance_matrix[col][row] += distance_matrix[row][col];
                distance_matrix[row][col] = 0;
            } else if (row == col) {
                distance_matrix[row][col] = 0;
            }
        }
    }
    return distance_matrix;
}

float *compute_divergence(float ** distance_matrix, std::vector<std::string> &names) {
    float divergence[names.size()];
    for (size_t num = 0; num < names.size(); num++) {
        divergence[num] = 0;
    }
    for (size_t row = 0; row < names.size(); row++) {
        for (size_t col = 0; col < names.size(); col++) {
            divergence[row] += distance_matrix[row][col];
            divergence[col] += distance_matrix[row][col];
        }
    }
    return divergence;
}

distance_matrix_grid make_new_distance_matrix_and_get_minimum(float **distance_matrix, float *divergence, std::vector<std::string> &names) {
    float **new_distance_matrix = new float * [names.size()];
    for (size_t i = 0; i < names.size(); i++) {
        new_distance_matrix[i] = new float [names.size()];
    }
    for (size_t row = 0; row < names.size(); row++) {
        for (size_t col = 0; col < names.size(); col++) {
            new_distance_matrix[row][col] = 0;
        }
    }
    for (size_t row = 0; row < names.size(); row++) {
        for (size_t col = 0; col < names.size(); col++) {
            if (row > col) {
                new_distance_matrix[row][col] = distance_matrix[row][col] - (divergence[row] + divergence[col]) / (names.size() - 2);
            }
        }
    }
    distance_matrix_grid grid = {1, 0};
    for (uint32_t row = 0; row < names.size(); row++) {
        for (uint32_t col = 0; col < names.size(); col++) {
            if (row > col && new_distance_matrix[grid.row][grid.col] > new_distance_matrix[row][col]) {
                grid = {row, col};
            }
        }
    }
    for (size_t i = 0; i < names.size(); i++) {
        delete [] new_distance_matrix[i];
    }
    delete [] new_distance_matrix;
    return grid;
}

float **distance_matrix_for_next_iteration(float **distance_matrix, float *divergence, distance_matrix_grid grid, std::vector<std::string> &names) {
    // get branch length for minimum grid
    float branch_length_row = distance_matrix[grid.row][grid.col] / 2 + (divergence[grid.row] - divergence[grid.col]) / (2 * (names.size() - 2));
    float branch_length_col = distance_matrix[grid.row][grid.col] - branch_length_row;
    std::string new_name = "(" + names[grid.row] + ":" + std::to_string(branch_length_row) + "," + names[grid.col] + ":" + std::to_string(branch_length_col) + ")";
    names.erase(names.begin() + grid.row);
    names.erase(names.begin() + grid.col);
    names.emplace_back(new_name);
    float **new_distance_matrix = new float * [names.size()];
    for (size_t i = 0; i < names.size(); i++) {
        new_distance_matrix[i] = new float [names.size()];
    }
    uint32_t row_pre = 0;
    uint32_t col_pre = 0;
    for (size_t row = 0; row < names.size() - 1; row++) {
        if (row_pre == grid.row || row_pre == grid.col) {
            row_pre++;
            if (row_pre == grid.row || row_pre == grid.col) {
                row_pre++;
            }
        }
        for (size_t col = 0; col < names.size(); col++) {
            if (col_pre == grid.col || col_pre == grid.row) {
                col_pre++;
                if (col_pre == grid.col || col_pre == grid.row) {
                    col_pre++;
                }
            }
            if (row <= col) {
                new_distance_matrix[row][col] = 0;
                col_pre++;
                if (col_pre >= names.size() + 1) {
                    col_pre = 0;
                }
            } else {
                new_distance_matrix[row][col] = distance_matrix[row_pre][col_pre];
                col_pre++;
            }
        }
        row_pre++;
    }
    col_pre = 0;
    for (size_t last_line; last_line < names.size() - 1; last_line++) {
        if (col_pre == grid.col || col_pre == grid.row) {
            col_pre++;
        }
        new_distance_matrix[names.size() - 1][last_line] = distance_matrix[std::max(col_pre, grid.col)][std::min(col_pre, grid.col)] + distance_matrix[std::max(col_pre, grid.row)][std::min(col_pre, grid.row)]
                 - distance_matrix[grid.row][grid.col] / 2;
        col_pre++;
    }
    new_distance_matrix[names.size() - 1][names.size() - 1] = 0;
    for (size_t i = 0; i < names.size(); i++) {
        delete [] distance_matrix[i];
    }
    delete [] distance_matrix;
    return new_distance_matrix;
}

void final_iteration(float **distance_matrix, std::vector<std::string> &names) {
    std::string new_name = "(" + names[0] + ":" + std::to_string(distance_matrix[1][0]) + " center," + names[1] + ")";
    names.clear();
    names.emplace_back(new_name);
}

void get_nj_tree(std::vector<std::string> &names, const char *file_path, int species_num) {
    float **distance_matrix = from_maf_to_distance_matrix(names, file_path, species_num);
    float *divergence = compute_divergence(distance_matrix, names);
    distance_matrix_grid grid = make_new_distance_matrix_and_get_minimum(distance_matrix, divergence, names);
    float **next_matrix = distance_matrix_for_next_iteration(distance_matrix, divergence, grid, names);
    while (names.size() != 2) {
        divergence = compute_divergence(next_matrix, names);
        grid = make_new_distance_matrix_and_get_minimum(next_matrix, divergence, names);
        next_matrix = distance_matrix_for_next_iteration(next_matrix, divergence, grid, names);
    }
    final_iteration(next_matrix, names);
}