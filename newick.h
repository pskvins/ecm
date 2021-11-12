//
// Created by 박석환 on 2021/08/09.
//
#include <gsl/gsl_linalg.h>

#ifndef ECM_EXPECTATION_STEP_H
#define ECM_EXPECTATION_STEP_H


class expectation_step {

};

//next is going to node of parent
struct newick_graph {
    newick_graph *next;
    std::string parent;
    double branch_length;
    std::string species;
    std::vector<newick_graph*> previous;
    double felsenstein[64];
    double upper[64];
    std::vector<int> zero_entries;
    bool base[64];
    gsl_matrix *expon_matrix;
    gsl_matrix *expon_eigen;
    gsl_matrix_complex *expon_eigen_com;
    double up_projection[64];
    gsl_complex up_projection_com[64];
    double down_projection[64];
    gsl_complex down_projection_com[64];
    int order;

    newick_graph() {
        next = NULL;
        parent = "None";
        branch_length = 0.0;
        species = "None";
        previous = {};
        order = 0;
        for (int num = 0; num < 64; num++) {
            felsenstein[num] = 0.0;
            upper[num] = 0.0;
            up_projection[num] = 0.0;
            down_projection[num] = 0.0;
            base[num] = true;
            up_projection_com[num].dat[0] = 0.0;
            up_projection_com[num].dat[1] = 0.0;
            down_projection_com[num].dat[0] = 0.0;
            down_projection_com[num].dat[1] = 0.0;
        }
        expon_matrix = gsl_matrix_alloc(64, 64);
        expon_eigen = gsl_matrix_alloc(64, 64);
        expon_eigen_com = gsl_matrix_complex_alloc(64, 64);

    }

    newick_graph(newick_graph *next, std::string parent, float branch_length, std::string species, std::vector<newick_graph*> previous) :
            previous(previous),
            parent(parent),
            branch_length(branch_length),
            species(species),
            next(next) {}

    void insert_inbetween_end_n_target(newick_graph *end, newick_graph *nodes, size_t count, newick_graph &child) {
        child.next = &nodes[count];
        nodes[count].previous.emplace_back(&child);
        nodes[count].next = end;
        for (size_t pos = 0; pos < end->previous.size(); pos++) {
            if (end->previous[pos] == &child) {
                end->previous.erase(end->previous.begin() + pos);
            }
        }
        if (std::find(end->previous.begin(), end->previous.end(), &nodes[count]) == end->previous.end()) {
            end->previous.emplace_back(&nodes[count]);
        }
    }
};

struct newick_start {
    std::vector<newick_graph*> next;

    newick_start() {
        next = {};
    }

    void connect_start(newick_start *start, newick_graph *current, size_t pos, newick_graph *end) {
        start->next.emplace_back(&current[pos]);
        current[pos].next = end;
        end->previous.emplace_back(&current[pos]);
    }
};

struct newick_iterator {
    newick_graph *now;
    newick_graph *parent;
    newick_start start_point;

    newick_iterator(newick_start &start_point) : start_point(start_point) {
        now = start_point.next[0];
        parent = now->next;
    }

    newick_iterator(newick_graph *now, newick_graph *parent, newick_start &start_point) :
            now(now),
            parent(parent),
            start_point(start_point) {}


    void to_start(size_t pos) {
        now = start_point.next[pos];
        parent = now->next;
    }

    void to_parent(void) {
        now = parent;
        parent = parent->next;
    }
};

#endif //ECM_EXPECTATION_STEP_H
