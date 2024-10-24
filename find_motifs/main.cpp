#include <iostream>

class DNA_data {
public:
    std::string genome_seq;
    std::string motif;
    int edit_dist;
    int motif_len;

    DNA_data() : genome_seq(""), motif(""), edit_dist(0),  motif_len(0) {}

    void set_data() {
        std::cout << "Enter edit distance k: " << std::endl;
        std::cin >> edit_dist;

        std::cout << "Enter DNA motif: " << std::endl;
        std::cin >> motif;
        motif_len = motif.length();

        std::cout << "Enter genome sequence: " << std::endl;
        std::cin >> genome_seq;
    }

    void find_motifs() {
        int motif_len = motif.length();
        int genome_len = genome_seq.length();

        for (int i = 0; i <= genome_len - motif_len; i++) {
            std::string genome_substring = genome_seq.substr(i, motif_len);
            int mismatches = 0;

            for (int j=0; j <= motif_len; j++) {
                if (motif[j] != genome_substring[j]) {
                    mismatches++;
                }
            }
            
            if (mismatches <= edit_dist) {
                std::cout << "Motif: " << motif << std::endl << "contains: " << mismatches << " mismatches";
            }
        }
    }

};

int main(){
    DNA_data dna;
    dna.set_data();
    dna.find_motifs();

    return 0;
}
