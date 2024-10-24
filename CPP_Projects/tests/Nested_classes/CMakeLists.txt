# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Project name
project(NestedClasses)

# Add executable target with all source files
add_executable(NestedClasses Main.cxx Outer.cxx)

# Add more debug information, including line number in case of segfault
set(CMAKE_BUILD_TYPE "Debug")
target_compile_options(NestedClasses PRIVATE -fno-omit-frame-pointer)
target_link_options(NestedClasses PRIVATE -fno-omit-frame-pointer)