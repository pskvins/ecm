#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

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

static std::vector<double> compute_pi_coding(const char * file_path_fasta, const char * file_path_gtf) {
    // read fasta file and divide into vector of name and sequence
    std::ifstream fasta;
    fasta.open(file_path_fasta);
    std::vector<std::string> name;
    std::vector<std::string> sequence;
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

    // read gtf file and get CDS info
    struct CDS_info {
        std::string chrom;
        uint32_t start;
        uint32_t end;
        uint additional;
        CDS_info(std::string chrom, uint32_t start, uint32_t end, uint additional):
            chrom(chrom),
            start(start),
            end(end),
            additional(additional){}
    };
    std::ifstream gtf;
    gtf.open(file_path_gtf);
    std::vector<CDS_info> all_CDS_info;
    std::string temp_chrom;
    uint32_t temp_start;
    uint32_t temp_end;
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
            }
            else if (queue == 2) {
                if (gtf_line.substr(0, pos) == "CDS") {
                    include = 1;
                }
                else {
                    include = 0;
                }
            }
            else if (queue == 3) {
                temp_start = static_cast<uint32_t>(std::stoul(gtf_line.substr(0, pos))) - 1;
            }
            else if (queue == 4) {
                temp_end = static_cast<uint32_t>(std::stoul(gtf_line.substr(0, pos)));
            }
            else if (queue == 7) {
                temp_additional = static_cast<uint>(std::stoul(gtf_line.substr(0, pos)));
            }
            queue++;
            gtf_line.erase(0, pos + delimiter.length());
        }
        if (include) {
            all_CDS_info.emplace_back(temp_chrom, temp_start, temp_end, temp_additional);
        }
        getline(gtf, gtf_line);
    }

    // extract codons
    std::vector<std::string> codon_set;
    for (size_t num = 0; num < all_CDS_info.size(); num++) {
        int pos_chrom;
        for (size_t pos = 0; pos < name.size(); pos++) {
            if (name[pos] == all_CDS_info[num].chrom) {
                break;
            }
            else {
                pos_chrom++;
            }
        }
        if (pos_chrom >= name.size()) {
            printf("No chromosome from gtf found from fasta");
        }
        for (uint32_t start = all_CDS_info[num].start + all_CDS_info[num].additional; start + 2 < all_CDS_info[num].end; start += 3) {
            std::string codon = "";
            for (int i = 0; i < 3; i++) {
                codon += sequence[pos_chrom][start + i];
            }
            codon_set.emplace_back(codon);
            codon = "";
            for (int j = 2; j >= 0; j--) {
                codon += complement_base(sequence[pos_chrom][start + j]);
            }
            codon_set.emplace_back(codon);
        }
    }
    std::vector<std::string> codon_list = {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"};
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
    std::vector<double> coding_codon_freq = compute_pi_coding("/Users/sukhwanpark/Downloads/sacCer3.fa", "/Users/sukhwanpark/Downloads/sacCer3.ncbiRefSeq.gtf");
    return 0;
}
