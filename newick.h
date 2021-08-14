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
    std::string species;
    std::vector<newick_graph*> previous;

    newick_graph() {
        next = NULL;
        parent = "None";
        branch_length = 0.0;
        species = "None";
        previous = {};
    }

    newick_graph(newick_graph *next, std::string parent, float branch_length, std::string species, std::vector<newick_graph*> previous) :
            previous(previous),
            parent(parent),
            branch_length(branch_length),
            species(species),
            next(next) {}

    void insert_inbetween_end_n_target(newick_graph &end, newick_graph *nodes, size_t count, newick_graph &child) {
        child.next = nodes + count;
        nodes[count].previous.emplace_back(&child);
        nodes[count].next = &end;
        for (size_t pos = 0; pos < end.previous.size(); pos++) {
            if (end.previous[pos] == &child) {
                end.previous.erase(end.previous.begin() + pos);
                break;
            }
        }
        end.previous.emplace_back(nodes + count);
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
    std::vector<newick_graph*> next;

    newick_start() {
        next = {};
    }

    newick_start(std::vector<newick_graph*> next) :
            next(next) {}

    void connect_start(newick_start &start, newick_graph *current, size_t pos, newick_graph &end) {
        start.next.emplace_back(current + pos);
        current[pos].next = &end;
        end.previous.emplace_back(current + pos);
    }
};

#endif //ECM_EXPECTATION_STEP_H
