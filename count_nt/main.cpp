#include <iostream>
#include <map>

class nt_counter {
private:
    std::map<char, int> nt_count_map;

public:
    // constructor to initialize the count
    nt_counter() {
        nt_count_map['A'] = 0;
        nt_count_map['C'] = 0;
        nt_count_map['G'] = 0;
        nt_count_map['T'] = 0;
    }

    void find_nt(const std::string& dna) {
        for (char nucleotide : dna) {
            if (nt_count_map.find(nucleotide) != nt_count_map.end()) {
                nt_count_map[nucleotide]++;
            }
        }
    }

    void print() {
        std::cout << "Nucleotide counts:\n";
        std::cout << "A: " << nt_count_map.at('A') << "\n";
        std::cout << "C: " << nt_count_map.at('C') << "\n";
        std::cout << "G: " << nt_count_map.at('G') << "\n";
        std::cout << "T: " << nt_count_map.at('T') << "\n";
    }
};

int main() {
    std::string dna = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC";

    nt_counter counter;
    counter.find_nt(dna);
    counter.print();

    return 0;
}
