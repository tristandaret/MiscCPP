# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Project name
project(MVEpavetext)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Tree)

# Include ROOT directories
include(${ROOT_USE_FILE})

# Add executable target
add_executable(MVEpavetext Main.cxx)

# Link ROOT libraries to the executable
target_link_libraries(MVEpavetext PUBLIC ${ROOT_LIBRARIES})

# Add more debug information, including line numbers in case of segfault
set(CMAKE_BUILD_TYPE "Debug")
target_compile_options(MVEpavetext PRIVATE -fno-omit-frame-pointer)
target_link_options(MVEpavetext PRIVATE -fno-omit-frame-pointer)
