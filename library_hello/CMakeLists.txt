cmake_minimum_required(VERSION 3.10)
project(HelloLibrary)

set(CMAKE_CXX_STANDARD 17)

# Define the shared library
add_library(HelloLibrary SHARED library.cpp)

# Add an executable for testing the library
add_executable(main_app main.cpp)

# Link the executable with the shared library
target_link_libraries(main_app PRIVATE HelloLibrary)



