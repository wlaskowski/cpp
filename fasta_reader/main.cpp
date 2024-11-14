#include <iostream>
#include <fstream>
#include <string>

class Sequence {
private:
    std::string sequence;

public:
    // initializing empty seq
    Sequence() : sequence("") {}

    // read a FASTA file
    void ReadFasta(const std::string &path) {
        std::ifstream file(path);
        if (!file.is_open()) {
            std::cerr << "Unable to open file: " << path << std::endl;
            return;
        }

        std::string line;
        while (std::getline(file, line)) {
            if (line.empty() or line[0] == '>') {
                continue;
            }
            sequence += line;
        }

        file.close();
    }

    // Method to print the stored sequence
    void PrintSeq() const {
        std::cout << "Sequence: " << std::endl << sequence << std::endl;
    }
};

int main() {
    Sequence seq;

    // Adjust the file path as needed for your system
    seq.ReadFasta("C:/Users/Wojtek/C++/untitled1/fasta.txt");

    // Print the loaded sequence
    seq.PrintSeq();

    return 0;
}
