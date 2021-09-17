//
// Created by 박석환 on 2021/08/15.
//

#include <iostream>
#include <vector>
#include <algorithm>

static std::string codon_list[64] = {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"};
double PhyloCSF_noncoding[64] = {0.050652, 0.017255, 0.021160, 0.032474, 0.018026, 0.008973, 0.008527, 0.015506, 0.018850, 0.010325, 0.010458, 0.015506, 0.035649, 0.014782, 0.017158, 0.032474, 0.019614, 0.010055, 0.009977, 0.017158, 0.010003, 0.006173, 0.005336, 0.010458, 0.008072, 0.006082, 0.005336, 0.008527, 0.014070, 0.009981, 0.009977, 0.021160, 0.022403, 0.007536, 0.009981, 0.014782, 0.011453, 0.006492, 0.006082, 0.010325, 0.010352, 0.006492, 0.006173, 0.008973, 0.016223, 0.007536, 0.010055, 0.017255, 0.028967, 0.016223, 0.014070, 0.035649, 0.017417, 0.010352, 0.008072, 0.018850, 0.017417, 0.011453, 0.010003, 0.018026, 0.028967, 0.022403, 0.019614, 0.050652};
double PhyloCSF_coding[64] = {0.042558, 0.024513, 0.030138, 0.036293, 0.017995, 0.012479, 0.008171, 0.020108, 0.020893, 0.010055, 0.009439, 0.014649, 0.018597, 0.016947, 0.020770, 0.030135, 0.026648, 0.007776, 0.012271, 0.013887, 0.017633, 0.006954, 0.005409, 0.013514, 0.003213, 0.002707, 0.001910, 0.006327, 0.013506, 0.005776, 0.010752, 0.012711, 0.045143, 0.020107, 0.019335, 0.037538, 0.016269, 0.012098, 0.006217, 0.020106, 0.011230, 0.009807, 0.006088, 0.022349, 0.012260, 0.011228, 0.010804, 0.021456, 0.001050, 0.014490, 0.000519, 0.019166, 0.019118, 0.014150, 0.008827, 0.023470, 0.000684, 0.005011, 0.010397, 0.008147, 0.026476, 0.018304, 0.026556, 0.026866};

struct CDS_info {
    std::string chrom;
    uint32_t start;
    uint32_t end;
    std::string strand;
    uint8_t additional;

    CDS_info(std::string chrom, uint32_t start, uint32_t end, std::string strand, uint8_t additional) :
            chrom(chrom),
            start(start),
            end(end),
            strand(strand),
            additional(additional) {}
};

std::string complement_base(char base) {
    if (base == 'A') {
        return "T";
    }
    else if (base == 'T') {
        return "A";
    }
    else if (base == 'C') {
        return "G";
    }
    else if (base == 'G') {
        return "C";
    }
    else if (base == '-') {
        return "-";
    }
    else if (base == 'N') {
        return "N";
    }
}

void return_complement_codon(std::string &codon) {
    std::string complement_codon = "";
    for (int num = 2; num >=0; num--) {
        complement_codon += complement_base(codon[num]);
    }
    codon = complement_codon;
}

void get_chr_and_seq(const char * file_path_fasta, std::vector<std::string>&name, std::vector<std::string>&sequence) {
    // read fasta file and divide into vector of name and sequence
    std::ifstream fasta;
    fasta.open(file_path_fasta);
    std::string line;
    getline(fasta, line);
    std::string delimiter = " ";
    size_t pos;
    while(!line.empty()) {
        if (line[0] == '>') {
            pos = line.find(delimiter);
            name.emplace_back(line.substr(1, pos - 1));
            getline(fasta, line);
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            sequence.emplace_back(line);
            getline(fasta, line);
            while (!line.empty() && line[0] != '>') {
                std::transform(line.begin(), line.end(), line.begin(), ::toupper);
                sequence.back() += line;
                getline(fasta, line);
            }
        }
    }
    fasta.close();
}

bool sort_CDS_info_vector(const CDS_info &first, const CDS_info &second) {
    if (first.chrom != second.chrom) {
        return (first.chrom < second.chrom);
    } else {
        return (first.start < second.start);
    }
}

std::vector<CDS_info> get_CDS_info(const char * file_path_gtf) {
    // read gtf file and get CDS info
    std::ifstream gtf;
    gtf.open(file_path_gtf);
    std::vector<CDS_info> all_CDS_info;
    std::string temp_chrom;
    uint32_t temp_start;
    uint32_t temp_end;
    std::string temp_strand;
    uint8_t temp_additional;
    std::string gtf_line;
    size_t pos;
    std::string delimiter = "\t";
    bool include = false;
    getline(gtf, gtf_line);
    if (gtf_line[0] == '#') {
        while (gtf_line[0] == '#') {
            getline(gtf, gtf_line);
        }
    }
    while (!gtf_line.empty()) {
        if (gtf_line[0] == '#') {
            while (gtf_line[0] == '#') {
                getline(gtf, gtf_line);
            }
            continue;
        }
        uint8_t queue = 0;
        while (queue < 8) {
            pos = gtf_line.find(delimiter);
            if (queue == 0) {
                temp_chrom = gtf_line.substr(0, pos);
            } else if (queue == 2) {
                if (gtf_line.substr(0, pos) == "CDS") {
                    include = true;
                } else {
                    include = false;
                }
            } else if (queue == 3) {
                temp_start = static_cast<uint32_t>(std::stoul(gtf_line.substr(0, pos))) - 1;
            } else if (queue == 4) {
                temp_end = static_cast<uint32_t>(std::stoul(gtf_line.substr(0, pos))) - 1;
            } else if (queue == 6) {
                temp_strand = gtf_line.substr(0, pos);
            } else if (queue == 7) {
                if (include) {
                    temp_additional = static_cast<uint8_t>(std::stoul(gtf_line.substr(0, pos)));
                }
            }
            queue++;
            gtf_line.erase(0, pos + delimiter.length());
        }
        if (include) {
            all_CDS_info.emplace_back(temp_chrom, temp_start, temp_end, temp_strand, temp_additional);
        }
        getline(gtf, gtf_line);
    }
    gtf.close();
    std::sort(all_CDS_info.begin(), all_CDS_info.end(), sort_CDS_info_vector);
    return all_CDS_info;
}

void get_codon(std::vector<std::string> &codon_set, std::string &codon, std::vector<std::string> &sequence, uint32_t pos_chrom, size_t pos) {
    codon = "";
    for (size_t i = 0; i < 3; i++) {
        codon += sequence[pos_chrom][pos + i];
    }
    for (int j = 0; j < 3; j++) {
        if (codon[j] == 'N' || codon[j] == '-') {
            codon = "";
        }
    }
    if (codon.size() > 0) {
        codon_set.emplace_back(codon);
    }
}

void get_complement_codon(std::vector<std::string> &codon_set, std::string &codon, std::vector<std::string> &sequence, uint32_t pos_chrom, size_t pos) {
    codon = "";
    for (int i = 2; i >= 0; i--) {
        codon += complement_base(sequence[pos_chrom][pos + i]);
    }
    for (int j = 0; j < 3; j++) {
        if (codon[j] == 'N' || codon[j] == '-') {
            codon = "";
        }
    }
    if (codon.size() > 0) {
        codon_set.emplace_back(codon);
    }
}

uint32_t finding_index(std::vector<std::string> &name, std::string chrom) {
    for (int num = 0; num < name.size(); num++) {
        if (name[num] == chrom) {
            return num;
        }
    }
    throw ("index not found");
}

std::vector<double> coding_codon_freq(std::vector<std::string> &name, std::vector<std::string> &sequence, std::vector<CDS_info> &all_CDS_info) {
    // extract codons and calculate codon frequencies
    std::vector<std::string> codon_set;
    std::string codon = "";
    for (size_t num = 0; num < all_CDS_info.size(); num++) {
        uint32_t pos_chrom = finding_index(name, all_CDS_info[num].chrom);
        if (pos_chrom >= name.size()) {
            printf("No chromosome from gtf found from fasta");
        }
        if (all_CDS_info[num].strand == "+") {
            for (uint32_t start = all_CDS_info[num].start + all_CDS_info[num].additional; start + 2 < all_CDS_info[num].end; start += 3) {
                for (int i = 0; i < 3; i++) {
                    codon += sequence[pos_chrom][start + i];
                }
                codon_set.emplace_back(codon);
                codon = "";
            }
        }
        else {
            for (uint32_t start = all_CDS_info[num].start + all_CDS_info[num].additional; start + 2 < all_CDS_info[num].end; start += 3) {
                std::string codon = "";
                for (int j = 2; j >= 0; j--) {
                    codon += complement_base(sequence[pos_chrom][start + j]);
                }
                codon_set.emplace_back(codon);
            }
        }
    }
    std::vector<double> codon_count;
    for (size_t codon_index = 0; codon_index < sizeof(codon_list) / sizeof(codon_list[0]); codon_index++) {
        codon_count.emplace_back(std::count(codon_set.begin(), codon_set.end(), codon_list[codon_index]));
    }
    double  codon_sum;
    for (size_t i = 0; i < codon_count.size(); i++) {
        codon_sum += codon_count[i];
    }
    for (size_t j = 0; j < codon_count.size(); j++) {
        codon_count[j] = codon_count[j] / codon_sum;
    }
    return codon_count;
}

bool compare_func(CDS_info first, CDS_info second) {
    return (first.chrom < second.chrom);
}

void mark_chrom_bool(std::vector<uint8_t> &chrom_bool, CDS_info elm_CDS_info) {
    for (uint32_t start = elm_CDS_info.start + elm_CDS_info.additional; start < elm_CDS_info.end; start++) {
        chrom_bool[start] = 1;
    }
}


std::vector<double> noncoding_codon_freq(std::vector<std::string> &name, std::vector<std::string> &sequence, std::vector<CDS_info> &all_CDS_info) {
    std::vector<std::string> codon_set;
    std::string codon = "";
    std::vector<std::string> new_sequences = {};
    std::string new_sequence = "";
    std::sort(all_CDS_info.begin(), all_CDS_info.end(), compare_func);
    uint32_t pos_chrom = finding_index(name, all_CDS_info[0].chrom);
    std::vector<uint8_t> chrom_bool(sequence[pos_chrom].size(), 0);
    for (size_t index = 0; index < all_CDS_info.size(); index++) {
        if ((index == 0 || all_CDS_info[index - 1].chrom == all_CDS_info[index].chrom) && index != all_CDS_info.size() - 1) {
            mark_chrom_bool(chrom_bool, all_CDS_info[index]);
        }
        else if (all_CDS_info[index - 1].chrom != all_CDS_info[index].chrom) {
            new_sequence = "";
            for (size_t pos = 0; pos < chrom_bool.size() - 2; pos++) {
                if (chrom_bool[pos] == 0 && chrom_bool[pos] != 'N' && chrom_bool[pos] != '-') {
                    new_sequence += sequence[pos_chrom][pos];
                }
            }
            new_sequences.emplace_back(new_sequence);
            pos_chrom = finding_index(name, all_CDS_info[index].chrom);
            chrom_bool.clear();
            chrom_bool.assign(sequence[pos_chrom].size(), 0);
            mark_chrom_bool(chrom_bool, all_CDS_info[index]);
        }
        else if (index == all_CDS_info.size() - 1) {
            new_sequence = "";
            for (size_t pos = 0; pos < chrom_bool.size() - 2; pos++) {
                if (chrom_bool[pos] == 0 && chrom_bool[pos] != 'N' && chrom_bool[pos] != '-') {
                    new_sequence += sequence[pos_chrom][pos];
                }
            }
            new_sequences.emplace_back(new_sequence);
        }
    }
    for (size_t big_size = 0; big_size < new_sequences.size(); big_size++) {
        for (size_t small_size = 0; small_size < new_sequences[big_size].size() - 2; small_size += 3) {
            get_codon(codon_set, codon, new_sequences, big_size, small_size);
            get_complement_codon(codon_set, codon, new_sequences, big_size, small_size);
        }
    }
    std::vector<double> codon_count;
    for (size_t codon_index = 0; codon_index < sizeof(codon_list) / sizeof(codon_list[0]); codon_index++) {
        codon_count.emplace_back(std::count(codon_set.begin(), codon_set.end(), codon_list[codon_index]));
    }
    double  codon_sum;
    for (size_t i = 0; i < codon_count.size(); i++) {
        codon_sum += codon_count[i];
    }
    for (size_t j = 0; j < codon_count.size(); j++) {
        codon_count[j] = codon_count[j] / codon_sum;
    }
    return codon_count;
}