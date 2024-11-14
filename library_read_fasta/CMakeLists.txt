cmake_minimum_required(VERSION 3.10)
project(read_fasta_lib)

set(CMAKE_CXX_STANDARD 17)

add_library(FastaRead SHARED fasta_read.cpp)

add_executable(read_fasta_lib main.cpp)

target_link_libraries(read_fasta_lib FastaRead)
