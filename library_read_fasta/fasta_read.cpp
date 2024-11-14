#include "fasta_read.h"
#include <fstream>
#include <iostream>

std::string FastaReader::readfasta(const std::string& filepath) {
  std::ifstream file(filepath);
  if (!file.is_open()) {
    std::cerr << "Unable to open file" << std::endl;
    return "";
  }

  std::string line;
  std::string sequence;

  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '>') {
      continue;
    }
    sequence += line;
  }

  file.close();
  return sequence;
}

