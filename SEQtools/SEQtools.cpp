#include "SEQtools.h"

// analyzer
std::map<char, int> Analyzer::ntCounter(const std::string& dna) {
    std::map<char, int> nucleotideCount = {{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};
    for (char nt : dna) {
        if (nucleotideCount.find(nt) != nucleotideCount.end()) {
            nucleotideCount[nt]++;
        }
    }
    return nucleotideCount;
}

std::string Analyzer::transcribeDNA(const std::string& dna) {
    std::string rna = dna;
    for (char& nt : rna) {
        if (nt == 'T') {
            nt = 'U';
        }
    }
    return rna;
}

std::string Analyzer::computeRevComp(const std::string& dna) {
    std::string reversed_complement = dna;
    std::reverse(reversed_complement.begin(), reversed_complement.end());

    for (char& nt : reversed_complement) {
        switch (nt) {
            case 'A': nt = 'T'; break;
            case 'T': nt = 'A'; break;
            case 'C': nt = 'G'; break;
            case 'G': nt = 'C'; break;
        }
    }
    return reversed_complement;
}

void Analyzer::printntCounter(const std::map<char, int>& counts) {
    std::cout << "Nucleotide counts: " << std::endl
              <<  "A: " << counts.at('A') << std::endl
              << "C: " << counts.at('C') << std::endl
              << "G: " << counts.at('G') << std::endl
              << "T: " << counts.at('T') << std::endl << "\n";
}

void Analyzer::printRNA(const std::string& rna) {
    std::cout << "RNA sequence: " << rna << std::endl << "\n";
}

void Analyzer::printcomputeRevComp(const std::string& revComp) {
    std::cout << "Reverse complement: " << revComp << std::endl << "\n";
}

// hamming
int Hamming::computeDistance(const std::string& s1, const std::string& s2) {
    if (s1.length() != s2.length()) {
        std::cerr << "hamming: DNA sequences are not equal" << std::endl << "\n";
        return -1;
    }

    int distance = 0;
    for (int i = 0; i < s1.length(); i++) {
        if (s1[i] != s2[i]) {
            distance++;
        }
    }
    return distance;
}

void Hamming::printDistance(int distance) {
    if (distance != -1) {
        std::cout << "Hamming distance: " << distance << std::endl << "\n";
    }
}

// obliczanie transition/transversion ratio
bool Hamming::isTransition(char a, char b) {
    return (a == 'A' && b == 'G') || (a == 'G' && b == 'A') ||
           (a == 'C' && b == 'T') || (a == 'T' && b == 'C');
}

double Hamming::computeRatio(const std::string& s1, const std::string& s2) {
    if (s1.length() != s2.length()) {
        std::cerr << "hamming: DNA sequences are not equal" << std::endl << "\n";
        return -1.0;
    }

    int transitions = 0, transversions = 0;

    for (int i = 0; i < s1.length(); i++) {
        if (s1[i] != s2[i]) {
            if (isTransition(s1[i], s2[i])) {
                transitions++;
            } else {
                transversions++;
            }
        }
    }

    if (transversions == 0 || transitions == 0) {
        std::cerr << "one of the mutations does not appear, ratio undefined" << std::endl << "\n";
        return -1.0;
    }

    return static_cast<double> (transitions) / static_cast<double>(transversions);
}

// translator
std::string Translator::transcribeTranslate(const std::string& dna) {
    std::string rna = transcribeDNA(dna);
    return translateRNA(rna);
}

std::string Translator::translateRNA(const std::string& rna) {
    const std::map<std::string, char> codonTable = {
        {"UUU", 'F'}, {"UUC", 'F'}, {"UUA", 'L'}, {"UUG", 'L'},
        {"UCU", 'S'}, {"UCC", 'S'}, {"UCA", 'S'}, {"UCG", 'S'},
        {"UAU", 'Y'}, {"UAC", 'Y'}, {"UAA", '*'}, {"UAG", '*'},
        {"UGU", 'C'}, {"UGC", 'C'}, {"UGA", '*'}, {"UGG", 'W'},
        {"CUU", 'L'}, {"CUC", 'L'}, {"CUA", 'L'}, {"CUG", 'L'},
        {"CCU", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
        {"CAU", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
        {"CGU", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
        {"AUU", 'I'}, {"AUC", 'I'}, {"AUA", 'I'}, {"AUG", 'M'},
        {"ACU", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
        {"AAU", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
        {"AGU", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
        {"GUU", 'V'}, {"GUC", 'V'}, {"GUA", 'V'}, {"GUG", 'V'},
        {"GCU", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
        {"GAU", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
        {"GGU", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
    };

    std::string protein;
    for (int i = 0; i + 2 < rna.size(); i += 3) {
        std::string codon = rna.substr(i, 3);
        auto codon_aa = codonTable.find(codon);
        if (codon_aa != codonTable.end() && codon_aa->second != '*') {
            protein += codon_aa->second;
        }
    }
    return protein;
}

void Translator::printProtein(const std::string& protein) {
    std::cout << "Protein sequence: " << protein << std::endl << "\n";
}

// splicing
std::string SpliceTranslate::removeIntrons(const std::string& dna, const std::vector<std::string>& introns) {
    std::string exon = dna;
    for (const std::string& intron : introns) {
        int pos;
        while ((pos = exon.find(intron)) != std::string::npos) {
            exon.erase(pos, intron.length());
        }
    }
    return exon;
}

std::string SpliceTranslate::processDNA(const std::string& dna, const std::vector<std::string>& introns) {
    std::string exons_joined = removeIntrons(dna, introns);
    return transcribeTranslate(exons_joined);
}

// ProtMass
std::map<char, double> ProtMass::getMonoisotopicMassTable() {
    return {
        {'A', 71.03711}, {'C', 103.00919}, {'D', 115.02694}, {'E', 129.04259},
        {'F', 147.06841}, {'G', 57.02146}, {'H', 137.05891}, {'I', 113.08406},
        {'K', 128.09496}, {'L', 113.08406}, {'M', 131.04049}, {'N', 114.04293},
        {'P', 97.05276}, {'Q', 128.05858}, {'R', 156.10111}, {'S', 87.03203},
        {'T', 101.04768}, {'V', 99.06841}, {'W', 186.07931}, {'Y', 163.06333}
    };
}

double ProtMass::computeProteinMass(const std::string& protein) {
    double totalMass = 0.0;
    std::map<char, double> massTable = getMonoisotopicMassTable();

    for (char aminoAcid : protein) {
        totalMass += massTable[aminoAcid];
    }

    return totalMass;
}

double ProtMass::computeMassFromDNA(const std::string& dna) {
    std::string protein = transcribeTranslate(dna);
    return computeProteinMass(protein);
}

