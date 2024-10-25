#include <iostream>

class DNA_data {
public:
    std::string genome_seq;
    std::string motif;
    int edit_dist;
    int motif_len;
    int genome_len;

    DNA_data() : genome_seq(""), motif(""), edit_dist(0),  motif_len(0), genome_len(0) {}

    void set_data() {
        std::cout << "Enter edit distance k: " << std::endl;
        std::cin >> edit_dist;

        std::cout << "Enter DNA motif: " << std::endl;
        std::cin >> motif;
        motif_len = motif.length();

        std::cout << "Enter genome sequence: " << std::endl;
        std::cin >> genome_seq;
        genome_len = genome_seq.length();
    }

    int calculate_distance(const std::string &fragm_a, const std::string &fragm_b) {
        int dist = 0;
        for (int i = 0; i < fragm_a.length(); i++) {
            if (fragm_a[i] != fragm_b[i]) {
                dist++;
            }
        }
        return dist;
    }

    void find_motifs() {
        for (int i = 0; i <= genome_len - motif_len; i++) {
            std::string genome_substring = genome_seq.substr(i, motif_len);

            int mismatches = calculate_distance(motif, genome_substring);

            // If mismatches are within the allowed edit distance, output the starting index and motif length
            if (mismatches <= edit_dist) {
                std::cout << i << " " << motif_len << std::endl;
            }
        }
    }
};

int main() {
    DNA_data dna;
    dna.set_data();
    dna.find_motifs();

    return 0;
}
