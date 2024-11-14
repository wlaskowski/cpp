#include "fasta_read.h"
#include <iostream>

int main() {
    FastaReader reader;

    std::string filepath = "C:/Users/Wojtek/C++/read_fasta_lib/example.fasta";
    std::string sequence = reader.readfasta(filepath);

    std::cout << "Sequence:\n" << sequence << std::endl;

    return 0;
}
