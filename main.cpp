#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

std::vector<std::string> codon_list = {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"};
std::vector<double> PhyloCSF = {0.0357983654718, 0.0146634321466, 0.0209916853044, 0.0241608300459, 0.0196605104347, 0.0113786862615, 0.0029009138093, 0.0163015785126, 0.0225609755057, 0.0145991965297, 0.0181825729424, 0.0164781957079, 0.0186856678405, 0.0130405682814, 0.0179309987722, 0.0240968720195, 0.0181775551111, 0.0146957289426, 0.0210180199751, 0.0178777837456, 0.0179375335248, 0.014026410045, 0.00304930001685, 0.0181980546371, 0.00237235102112, 0.00247910724881, 0.00304553476624, 0.00292174613484, 0.0129216428963, 0.016886060241, 0.0212194293029, 0.021001239545, 0.0207025677211, 0.0101079314077, 0.0168739180256, 0.0130368421188, 0.0145312058205, 0.0117191223335, 0.00248094305692, 0.0145761834505, 0.0160538354021, 0.0117605522361, 0.0140380390381, 0.0114567437712, 0.0112780796487, 0.0101486396189, 0.0148519084459, 0.0147005486842, 0.0206773071958, 0.0111337585868, 0.012890738691, 0.0187059909187, 0.0197741432688, 0.0160489156219, 0.00237751747496, 0.0226108300737, 0.0197730187728, 0.0145713864549, 0.0179588059761, 0.0198394530312, 0.0208205335746, 0.0208303362555, 0.0183984258157, 0.0360132307673};

struct CDS_info {
    std::string chrom;
    uint32_t start;
    uint32_t end;
    std::string strand;
    uint additional;

    CDS_info(std::string chrom, uint32_t start, uint32_t end, std::string strand, uint additional) :
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

void get_chr_and_seq(const char * file_path_fasta, std::vector<std::string>&name, std::vector<std::string>&sequence) {
    // read fasta file and divide into vector of name and sequence
    std::ifstream fasta;
    fasta.open(file_path_fasta);
    std::string line;
    getline(fasta, line);
    while(line.size() > 0) {
        if (line[0] == '>') {
            name.emplace_back(line);
            getline(fasta, line);
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            sequence.emplace_back(line);
            getline(fasta, line);
            while (line.size() > 0 && line[0] != '>') {
                std::transform(line.begin(), line.end(), line.begin(), ::toupper);
                sequence.back() += line;
                getline(fasta, line);
            }
        }
    }
    fasta.close();
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
    uint temp_additional;
    std::string gtf_line;
    size_t pos = 0;
    std::string delimiter = "\t";
    bool include = 0;
    getline(gtf, gtf_line);
    while (gtf_line.size() > 0) {
        uint queue = 0;
        while (queue < 8) {
            pos = gtf_line.find(delimiter);
            if (queue == 0) {
                temp_chrom = gtf_line.substr(0, pos);
            } else if (queue == 2) {
                if (gtf_line.substr(0, pos) == "CDS") {
                    include = 1;
                } else {
                    include = 0;
                }
            } else if (queue == 3) {
                temp_start = static_cast<uint32_t>(std::stoul(gtf_line.substr(0, pos))) - 1;
            } else if (queue == 4) {
                temp_end = static_cast<uint32_t>(std::stoul(gtf_line.substr(0, pos)));
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
    return all_CDS_info;
}

uint32_t finding_index(std::vector<std::string> &name,std::string chrom) {
    int pos_chrom = 0;
    for (size_t pos = 0; pos < name.size(); pos++) {
        if (name[pos] == ">" + chrom) {
            break;
        }
        else {
            pos_chrom++;
        }
    }
    return pos_chrom;
}

std::vector<double> coding_codon_freq(std::vector<std::string> &name, std::vector<std::string> &sequence, std::vector<CDS_info> &all_CDS_info) {
    // extract codons and calculate codon frequencies
    std::vector<std::string> codon_set;
    for (size_t num = 0; num < all_CDS_info.size(); num++) {
        uint32_t pos_chrom = finding_index(name, all_CDS_info[num].chrom);
        if (pos_chrom >= name.size()) {
            printf("No chromosome from gtf found from fasta");
        }
        if (all_CDS_info[num].strand == "+") {
            for (uint32_t start = all_CDS_info[num].start + all_CDS_info[num].additional; start + 2 < all_CDS_info[num].end; start += 3) {
                std::string codon = "";
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
    for (size_t codon_index = 0; codon_index < codon_list.size(); codon_index++) {
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
    if (elm_CDS_info.strand == "+") {
        for (uint32_t start = elm_CDS_info.start; start < elm_CDS_info.end; start++) {
            if (chrom_bool[start] == 0) {
                chrom_bool[start] = 1;
            }
            else if (chrom_bool[start] == 2) {
                chrom_bool[start] =3;
            }
        }
    }
    else {
        for (uint32_t start = elm_CDS_info.start; start < elm_CDS_info.end; start++) {
            if (chrom_bool[start] == 0) {
                chrom_bool[start] = 2;
            }
            else if (chrom_bool[start] == 1) {
                chrom_bool[start] =3;
            }
        }
    }
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

std::vector<double> noncoding_codon_freq(std::vector<std::string> &name, std::vector<std::string> &sequence, std::vector<CDS_info> &all_CDS_info) {
    std::vector<std::string> codon_set;
    std::string codon = "";
    std::sort(all_CDS_info.begin(), all_CDS_info.end(), compare_func);
    uint32_t pos_chrom = finding_index(name, all_CDS_info[0].chrom);
    std::vector<uint8_t> chrom_bool(sequence[pos_chrom].size(), 0);
    for (size_t index = 0; index < all_CDS_info.size(); index++) {
        if (index == 0 || all_CDS_info[index - 1].chrom == all_CDS_info[index].chrom) {
            mark_chrom_bool(chrom_bool, all_CDS_info[index]);
        }
        else if (all_CDS_info[index - 1].chrom != all_CDS_info[index].chrom) {
            for (size_t pos = 0; pos < chrom_bool.size() - 2; pos++) {
                if (chrom_bool[pos] == 3 || chrom_bool[pos + 1] == 3 || chrom_bool[pos+2] == 3) {
                    while (chrom_bool[pos + 1] == 3 || chrom_bool[pos + 2] == 3) {
                        pos++;
                    }
                }
                else if (chrom_bool[pos] == 0 && chrom_bool[pos] == chrom_bool[pos + 1] && chrom_bool[pos] == chrom_bool[pos + 2]) {
                    get_codon(codon_set, codon, sequence, pos_chrom, pos);
                    get_complement_codon(codon_set, codon, sequence, pos_chrom, pos);
                }
                else if ((chrom_bool[pos] == 0 || chrom_bool[pos] == 2) || (chrom_bool[pos + 1] == 0 || chrom_bool[pos + 1] == 2) || (chrom_bool[pos + 2] == 0 || chrom_bool[pos + 2] == 2)) {
                    get_codon(codon_set, codon, sequence, pos_chrom, pos);
                }
                else if ((chrom_bool[pos] == 0 || chrom_bool[pos] == 1) || (chrom_bool[pos + 1] == 0 || chrom_bool[pos + 1] == 1) || (chrom_bool[pos + 2] == 0 || chrom_bool[pos + 2] == 1)) {
                    get_complement_codon(codon_set, codon, sequence, pos_chrom, pos);
                }
            }
            pos_chrom = finding_index(name, all_CDS_info[index].chrom);
            chrom_bool.clear();
            chrom_bool.assign(sequence[pos_chrom].size(), 0);
            mark_chrom_bool(chrom_bool, all_CDS_info[index]);
        }
    }
    std::vector<double> codon_count;
    for (size_t codon_index = 0; codon_index < codon_list.size(); codon_index++) {
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


int main() {
    std::vector<std::string> name;
    std::vector<std::string> sequence;

    get_chr_and_seq("/Users/sukhwanpark/Downloads/sacCer3.fa", name, sequence);

    std::vector<CDS_info> all_CDS_info = get_CDS_info("/Users/sukhwanpark/Downloads/sacCer3.ncbiRefSeq.gtf");

    std::vector<double> noncoding_freq = noncoding_codon_freq(name, sequence, all_CDS_info);

    std::vector<double> diff;

    for (size_t i = 0; i < noncoding_freq.size(); i++) {
        diff.emplace_back(PhyloCSF[i] - noncoding_freq[i]);
    }

    return 0;
}
