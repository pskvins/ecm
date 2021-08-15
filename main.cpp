#include <iostream>
#include <vector>
#include <algorithm>
#include "njtree.cpp"
#include "newick.cpp"
#include "codon_frequency.cpp"


int main() {
    newick_start start_point;
    newick_graph end_point;

    newick_graph *nodes = process_newick("/Users/sukhwanpark/Downloads/7yeast.nh", start_point, end_point);




    return 0;
}
