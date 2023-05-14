cmake_minimum_required(VERSION 3.16)

project(cnn)

include_directories(/opt/homebrew/Cellar/openblas/0.3.23/include)

link_libraries("/opt/homebrew/Cellar/openblas/0.3.23/lib/libopenblas.dylib")

add_executable(cnn test.cpp)