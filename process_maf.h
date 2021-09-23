//
// Created by 박석환 on 2021/08/15.
//

#ifndef ECM_PROCESS_MAF_H
#define ECM_PROCESS_MAF_H


class process_maf {

};

int codonTable[4][4][4] = {{{0, 1, 3, 2}, {4, 5, 7, 6}, {8, 9, 11, 10}, {12, 13, 15,14}},
                           {{16, 17, 19, 18}, {20, 21, 23, 22}, {24, 25, 27, 26}, {28, 29, 31, 30}},
                           {{32, 33, 35, 34}, {36, 37, 39, 38}, {40, 41, 43, 42}, {44, 45, 47, 46}},
                           {{48, 49, 51, 50}, {52, 53, 55, 54}, {56, 57, 59, 58}, {60, 61, 63, 62}}};

struct aligned_codon {
    std::string species;
    char codon;

    aligned_codon() {
        species = "None";
        codon = 64;
    }

    aligned_codon(std::string species, char codon) :
        species(species),
        codon(codon) {}
};

#endif //ECM_PROCESS_MAF_H
