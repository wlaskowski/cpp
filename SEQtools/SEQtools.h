#ifndef SEQTOOLS_H
#define SEQTOOLS_H

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <algorithm>


// nt counter, DNA to RNA, DNA rev comp
class Analyzer {
public:
    std::map<char, int> ntCounter(const std::string& dna);
    std::string transcribeDNA(const std::string& dna);
    std::string computeRevComp(const std::string& dna);

    void printntCounter(const std::map<char, int>& counts);
    void printRNA(const std::string& rna);
    void printcomputeRevComp(const std::string& revComp);
};

// HAMM, TRAN
class Hamming {
public:
    int computeDistance(const std::string& s, const std::string& t);
    void printDistance(int distance);

    double computeRatio(const std::string& s, const std::string& t);
private:
    bool isTransition(char a, char b);
};

// DNA to protein
class Translator : public Analyzer {
public:
    std::string transcribeTranslate(const std::string& dna);
    void printProtein(const std::string& protein);

protected:
    std::string translateRNA(const std::string& rna);
};

// splicing
class SpliceTranslate : public Translator {
public:
    std::string removeIntrons(const std::string& dna, const std::vector<std::string>& introns);
    std::string processDNA(const std::string& dna, const std::vector<std::string>& introns);
};

// calculating prot mass
class ProtMass : public Translator {
public:
    double computeProteinMass(const std::string& protein);
    double computeMassFromDNA(const std::string& dna);
protected:
    std::map<char, double> getMonoisotopicMassTable();
};


#endif // SEQTOOLS_H
