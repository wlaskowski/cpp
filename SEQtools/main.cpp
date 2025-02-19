#include <iostream>
#include "SEQtools.h"

int main() {
    Analyzer analyzer; // nt counting, DNA to RNA, rev comp DNA
    Hamming hamming; // hamming + calculating transition/transversion ratio
    Translator translator; // DNA to RNA to protein
    SpliceTranslate spliceTranslate; // DNA -> splicing -> RNA to protein
    ProtMass protMass; // calculating protein mass

    // DNA sequence
    std::string dna = "ATGTCCTGTGATGCAAGTACCGTGTTATCCGGGCTATCGTCGAATGGCCCTCTTAGTTGGTACGCTCTTTATCATTCCATTGAGACATCGTCAGCACTTCACAGCGTCAACTATCCCGAAAGGCGGCCCCCAAACCTTTTATTTAATACATAAGGCGCTCTAATTCCAAGCACGAAAAATCATAGGGCCTTCGACAGACGCCGCATATTGCTGCCGTTAAATATTACCATGTACACTTCACCCGGCTACAGCTTGTAGCACAGGAAATTTTACCGCCGACGGAGAGCTTGAGACGCGTGAGTTTAACGTCTGGATGTGGCAGCCGGTGCCAGGGCCGGTTTGCAGTACTGGTACTTGGCATGCCTCGGGCGTCATGCGTAAGTCACGACCGCTTGCCAACCAGTATCCTTATGAGCTTGTGTATCGCTATTGGGGAACCTCTCCCACATGTTCGGAAGTCCAAGCAAACGGTATGGGTCATGATTGCTCTAACGGCCCTGAACATAGTAATCAGTACGGAAACTCTTACGTCCTCTTGAGCATATGATGGTCGATTATCGACTTAGCTGAGTGCTGGTATTAACTGTCTGTCGTCCCATGAGTTTCAGAAGGGAAGTTGGAGGAGCTTAAGACTACCCCCTAGGCTAATATTGCGGTCCTCCCGTAGCAATGTACGGCCATCGCTGTGAATTCAGAGAAAAGAGGGCTTACTCTTTGGCTTCACCGTTGAGGGTGAAACCCTTACCCGTCGAGGACCCCCCTATACGGCCGGTGCGCGCGGGCGTGGCTCCGGGAGTCCGACCGAATGTCAGCGGTTAGGCCATATAGTCAGTATAATTCTTGACACTTTCTATTTTTACAGTACTCCAAATCGTTACATTCTTTGCTTGAAGATCAGTCTGTGCAATTTTATGACTGGCACCTTCGACTCTACTAG";

    // intron sequences
    std::vector<std::string> introns = {
        "AAAAGAGGGCTTACTCTTTGGCT",
    "GAAATTTTACCGCCGACGGAGAGCTTG",
    "TGAAGATCAGTCTGTGCAA",
    "TTGCCAACCAGTATCCTTATGAGCTTGTGTATCGCTATTGGGGAA",
    "CGGCCGGTGCGCGCGGGCGTGGCTCCGGGAGTCCGACCGAATG",
    "TATAATTCTTGACACTTTCTATTTTTACAGTACTCCAAATC",
    "TATTGCTGCCGTTAAATATTACCATGTACACTTCACCCGGCTA",
    "TTAAGACTACCCCCTAGGCTAATATTGCGGTCCTCCCGTAG",
    "TGCTGGTATTAACTGTCTGTCGTCCCATGAGTTTCAG",
    "CAAACGGTATGGGTCATGATTGCTCTAACGGCCCTGAACATAGTA",
    "GGATGTGGCAGCCGGTGCCAG",
    "GCTCTTTATCATTCCATTGAGACATCGT",
    "TGAGCATATGATGG",
    "ACCTTTTATTTAATACATAAGGCGCTCTAA"
    };

    // Hamming distance and mutation ratio sequences
    std::string dna1 = "GGAGAGCCCAGGTCTGGCAGGCCAAACATCCATTGCAGACGCAAAAGTAAGCTCTATTCTCTCGTGACGATCGGGCGGCGACATCTAGATCTTAGCTAGCAGGCAGAGTATCCTATTAGCTTAACCGTCCTCGAAGTCGAGCATATAGAACAGTTGTGATGCACGTTTCGCTATATTTAACATCGGTACTTAGGCCCGCCACCGTCCAGCCGTGATCGTTAAGGATGTAAATGGCCTCGGTCAATTAAAGAGAACACAGCACGGCCGCAAAGTTGTTGGGGGCCCGTTCGGTACACGAGGTACGACCATCCACCTAGTACGAGTTCGGTTCCAACATTTCGCTCCACCCTGGTTAAGCCGTTCGTTTTGTACCATGGCAGCGGACTAACGGTACCTGAAGAGGCTGAAAGAGTGCTGCCTACTGCCTACGGAAATTACGTTATGGCATACTCCCACGTCCACTGCACGCCATCTTCGTTTAAGGTCGGCAAAATGAGTTGTTGCTGAGCAGATAAGAAGGGCGTTAGATTGGATGGTCGCGGGTGACATTCGTCAACCTAAGGAATTGATAAAAGAATTTCGAACCTCATTACCAACCCCGCTTTCCTGATGCAGCCTGATTCAAGTGGACGAATACTGACGAACTCTCGGCATTCAAAACCAGACAGTGATTCCTTTGCTTCTTACGCTCAAATCAGCACGTGCCCTGGGGCATTCGCACCTGCCCCATCCTCGGCACTACACAGGATGCCCCAAATTGAGTTTTGTCAGGAAGATTCTGTTAGCACTATAATCTTTTATTGCAACCCTGAGCACACTGTAGAAGTGTATACCCTCATACGTCGGAAACAGATGTACGAAATTAGGGGACTAACAAATAGCGCTGCATGGATCTCAGTTGCAATGGCTCCGCGTAAGAGGATGTTAATACGCACAGGGACATACCGCGAGTC";
    std::string dna2 = "CCAGAGCTCGGGTTGGATAGGCCAGATTTCTACTGCAAGCGCCGAAGTAAGGTCTGGTACCTTAGCATGAATTGGCAGCGCTACCCGGGTCGTGACTAGTAGGTAGAGTGCCCCGTCGGCTTAACAGTCTTTGAAGTCGCGCATGTAGAGTAGTTGTGCTGTACGTTTCACCACATTTACTATCGGCACTTAGGCCGGCCGCTGCGTAGCCATGATTGATAAGGCTGTAAACGGCCTCGGTTGACTAAATAGCATATAAGACGACCGCAGAGCTGTCGGGAGCTCCTTCTATGCATGGGGTGGAACAGTCCAGCTAGCGCAAGTTCGACTCTGACATTTCAGTTCACCCTGGTTAAGCCGTTCATTCTGTACCGTGACAGAGGACAAACGGTACCTGGCGAAGATCAGGCAATACCTCTTATGGCTCATGAAAATTGCGTTATCCCATACTACCCCGTACACTGCACGCTACTTTCCTGCGGGGCCGGCGAAATGCGCTATTGTTATGGCGTTGAGAAGGGCGTCACGTTGGGTGGGTGCCAGTGGCGGCCGAAAACCTGAAGAAGTGATAAGGGAATTTCGAACCTCAACGCCAGCCCCCCTTCTCTAGATCAGCATGATTCAAGTGCACGGATACTGACGGACTTTCAGCCCTCAAGACTAGACGATGAGTCCATTATCTCATGCGATCAAACAGGCACATCCCCCGGAGCCTACGCACCTGCTCCATCCTCGGCGCTACGCGGGGTAGCTTAAATTGAGTTTCGTTAGGAAGATCCCGTCAGTATTGCAATCTTTGATTGCAGTTCTGAGTTCACTGTAGAGGTATATACTCTCAGGTGTTAGAGGGTGATCTATGAAGTTGGAGAACTAACGGGCCGCGCAGCATGGATCTTATTTGAAATGGTCCCGCGTGAAAAGGTGTTGTTACGCACAGAGACAGATTACGAGCC";

    // nt counting
    std::map<char, int> nt_counts = analyzer.ntCounter(dna);
    analyzer.printntCounter(nt_counts);

    // DNA to RNA
    std::string rna = analyzer.transcribeDNA(dna);
    analyzer.printRNA(rna);

    // rev comp DNA
    std::string revComp = analyzer.computeRevComp(dna);
    analyzer.printcomputeRevComp(revComp);

    // Hamming count
    int distance = hamming.computeDistance(dna1, dna2);
    hamming.printDistance(distance);

    // mutation ratio
    double ratio = hamming.computeRatio(dna1, dna2);
    if (ratio != -1.0) {
        std::cout << "Transition/Transversion Ratio: " << ratio << std::endl << "\n";
    }

    // DNA to protein
    std::string proteinPROT = translator.transcribeTranslate(dna);
    translator.printProtein(proteinPROT);

    // splicing
    std::string proteinSPLC = spliceTranslate.processDNA(dna, introns);
    spliceTranslate.printProtein(proteinSPLC);

    double prot_mass_from_dna = protMass.computeMassFromDNA(dna);
    std::cout << "protein mass: " << prot_mass_from_dna << " Da" << std::endl;

    return 0;
}
