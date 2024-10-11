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
};

int main(){
    DNA_data dna;
    dna.set_data();

    return 0;
}
