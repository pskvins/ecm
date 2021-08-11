//
// Created by 박석환 on 2021/08/09.
//

#ifndef ECM_EXPECTATION_STEP_H
#define ECM_EXPECTATION_STEP_H


class expectation_step {

};

//next is going to node of parent
struct newick_graph {
    newick_graph *next;
    std::string parent;
    float branch_length;
    std::string child;
    uint16_t rank;
    std::vector<newick_graph*> previous;

    newick_graph(){
        next = NULL;
        parent = "Noone";
        branch_length = 0.0;
        child = "No";
        rank = 0;
    }

    newick_graph(newick_graph *next, std::string parent, float branch_length, std::string child, uint16_t rank, std::vector<newick_graph*> previous) :
            previous(previous),
            parent(parent),
            branch_length(branch_length),
            child(child),
            rank(rank),
            next(next) {}

    void set_next(newick_graph *node) {
        next = node;
    }

    void set_previous(newick_graph *node) {
        previous.emplace_back(node);
    }

    void set_parent(std::string name) {
        parent = name;
    }

    void set_bl(float length) {
        branch_length = length;
    }

    void connect_end(newick_graph end, newick_graph *current) {
        end.previous.emplace_back(current);
    }

    void insert_inbetween_end_n_target(newick_graph end, newick_graph current, newick_graph child) {
        child.next = &current;
        current.previous.emplace_back(&child);
        current.next = &end;
        for (size_t pos = 0; pos < end.previous.size(); pos++) {
            if (end.previous[pos] == &child) {
                end.previous.erase(end.previous.begin() + pos);
                break;
            }
        }
        end.previous.emplace_back(&current);
    }
};

struct newick_iterator {
    newick_graph *now;
    newick_graph *parent;

    newick_iterator(newick_graph *now, newick_graph *parent) :
            now(now),
            parent(parent) {}

    void to_parent(void) {
        now = parent;
        parent = parent->next;
    }

    void to_end(newick_graph *end) {
        parent = end;
        now = end;
    }
};

struct newick_start {
    std::string parent;
    std::vector<newick_graph*> next;

    newick_start() :
    parent(NULL) {}

    newick_start(std::string parent, std::vector<newick_graph*> next) :
            parent(parent),
            next(next) {}

    void connect_start(newick_graph *current) {
        next.emplace_back(current);
    }
};

#endif //ECM_EXPECTATION_STEP_H
