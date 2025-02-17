cmake_minimum_required(VERSION 3.29)
project(rosalindtree)

set(CMAKE_CXX_STANDARD 20)

add_executable(rosalindtree main.cpp
        SEQtools.h
        SEQtools.cpp
)
