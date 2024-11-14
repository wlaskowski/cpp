#ifndef FASTA_READ_H
#define FASTA_READ_H

#include <string>

class FastaReader {
public:
  std::string readfasta(const std::string& filepath);
  };
#endif
